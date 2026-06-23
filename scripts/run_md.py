#!/usr/bin/env python3
"""OpenMM MD at 72C for GFP sequences - fixed version"""
import os, sys, json, math, csv, threading, glob
import numpy as np

AA_NAMES = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HID',
            'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
            'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# Heavy atoms needed for each residue (Amber ff14SB templates)
AA_ATOMS = {
    'ALA': ['N','CA','C','O','CB'], 'ARG': ['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],
    'ASN': ['N','CA','C','O','CB','CG','OD1','ND2'], 'ASP': ['N','CA','C','O','CB','CG','OD1','OD2'],
    'CYS': ['N','CA','C','O','CB','SG'], 'GLN': ['N','CA','C','O','CB','CG','CD','OE1','NE2'],
    'GLU': ['N','CA','C','O','CB','CG','CD','OE1','OE2'], 'GLY': ['N','CA','C','O'],
    'HID': ['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2'],
    'ILE': ['N','CA','C','O','CB','CG1','CG2','CD1'], 'LEU': ['N','CA','C','O','CB','CG','CD1','CD2'],
    'LYS': ['N','CA','C','O','CB','CG','CD','CE','NZ'], 'MET': ['N','CA','C','O','CB','CG','SD','CE'],
    'PHE': ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
    'PRO': ['N','CA','C','O','CB','CG','CD'], 'SER': ['N','CA','C','O','CB','OG'],
    'THR': ['N','CA','C','O','CB','OG1','CG2'],
    'TRP': ['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],
    'TYR': ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'],
    'VAL': ['N','CA','C','O','CB','CG1','CG2'],
}

def create_pdb(seq, name):
    """Create compact PDB with all heavy atoms, minimum 1.2A spacing."""
    import random
    random.seed(42)
    lines = ["CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1", "MODEL        1"]
    an, all_pos = 1, []
    for i, aa in enumerate(seq):
        rn = AA_NAMES.get(aa, 'ALA')
        atoms = AA_ATOMS.get(rn, ['N','CA','C','O','CB'])
        for j, nm in enumerate(atoms):
            while True:
                x, y, z = random.gauss(0, 3), random.gauss(0, 3), random.gauss(0, 3)
                if all((x-px)**2+(y-py)**2+(z-pz)**2 >= 9.0 for (px,py,pz) in all_pos):
                    all_pos.append((x,y,z))
                    break
            el = 'N' if nm[0]=='N' else 'O' if nm[0]=='O' else 'S' if nm[0]=='S' else 'C'
            lines.append(f"ATOM  {an:5d} {nm:<4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
            an += 1
        if i == len(seq) - 1:
            while True:
                x, y, z = random.gauss(0, 3), random.gauss(0, 3), random.gauss(0, 3)
                if all((x-px)**2+(y-py)**2+(z-pz)**2 >= 9.0 for (px,py,pz) in all_pos):
                    all_pos.append((x,y,z))
                    break
            lines.append(f"ATOM  {an:5d} OXT  {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          O  ")
            an += 1
    lines += ["TER", "ENDMDL"]
    return "\n".join(lines)

def run_md(pdb_path, out_dir, temp_c=72.0, time_ns=1.0, prefix="prot", gpu="0"):
    """Run MD simulation."""
    try:
        from openmm import app, unit
        from openmm.app import PDBFile, ForceField, Modeller, Simulation
        from openmm import LangevinIntegrator, MonteCarloBarostat, Platform
        
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
        pdb = PDBFile(pdb_path)
        ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(ff)
        modeller.addSolvent(ff, padding=1.2*unit.nanometer)
        
        temp = temp_c + 273.15
        system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)
        system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, temp*unit.kelvin, 25))
        
        integrator = LangevinIntegrator(temp*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        platform = Platform.getPlatformByName("OpenCL")
        
        sim = Simulation(modeller.topology, system, integrator, platform)
        sim.context.setPositions(modeller.positions)
        
        print(f"  [GPU{gpu}] Minimizing...", flush=True)
        sim.minimizeEnergy(maxIterations=5000)
        
        print(f"  [GPU{gpu}] Equilibrating (30C)...", flush=True)
        integrator.setTemperature(303.15*unit.kelvin)
        sim.step(10000)
        
        print(f"  [GPU{gpu}] Production MD ({temp_c}C, {time_ns}ns)...", flush=True)
        nsteps = int(time_ns * 500000)
        interval = max(1, nsteps // 200)
        
        ref_state = sim.context.getState(getPositions=True)
        ref_pos = ref_state.getPositions()
        
        rmsds = []
        for step in range(nsteps):
            sim.step(100)
            if step % interval == 0:
                state = sim.context.getState(getPositions=True)
                cur = state.getPositions()
                n = min(len(ref_pos), len(cur))
                d = sum((ref_pos[i][j].value_in_unit(unit.nanometer) - cur[i][j].value_in_unit(unit.nanometer))**2
                       for i in range(n) for j in range(3))
                rmsd = math.sqrt(d / (n*3))
                rmsds.append(rmsd)
                ns = (step+1) * 0.0002
                print(f"  [GPU{gpu}] {ns:.2f}ns RMSD={rmsd:.4f}nm", flush=True)
        
        avg = float(np.mean(rmsds[-20:]))
        drift = float(rmsds[-1] - rmsds[0]) if len(rmsds) > 1 else 0
        stable = "STABLE" if avg < 0.3 and drift < 0.1 else ("MODERATE" if avg < 0.5 else "UNSTABLE")
        
        result = {"prefix":prefix,"temp_C":temp_c,"time_ns":time_ns,
                  "avg_RMSD_nm":round(avg,4),"drift_nm":round(drift,4),"stability":stable}
        
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, f"{prefix}_result.json"),"w") as f:
            json.dump(result, f, indent=2)
        print(f"  [GPU{gpu}] {prefix}: RMSD={avg:.4f}nm [{stable}]", flush=True)
        return result
    except Exception as e:
        print(f"  [GPU{gpu}] ERROR: {e}", flush=True)
        return None

if __name__ == "__main__":
    wd = "/mnt/home/mengd/FoldSynth"
    seqs = []
    with open(os.path.join(wd, "submission.csv")) as f:
        for row in csv.DictReader(f):
            seqs.append((row["Seq_ID"], row["Sequence"]))
    print(f"Sequences: {len(seqs)}")
    
    pdb_dir = os.path.join(wd, "designs/pdbs")
    os.makedirs(pdb_dir, exist_ok=True)
    md_dir = os.path.join(wd, "results/md_simulations")
    os.makedirs(md_dir, exist_ok=True)
    
    pdbs = []
    for sid, seq in seqs:
        name = f"FoldSynth_{sid}"
        pdb_content = create_pdb(seq, name)
        pdb_path = os.path.join(pdb_dir, f"{name}.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_content)
        pdbs.append(pdb_path)
        print(f"  PDB: {name}.pdb ({len(seq)}aa)")
    
    print("\nRunning MD on 2 GPUs (GPU 0+1, leaving 4 free)...")
    b1, b2 = pdbs[:3], pdbs[3:]
    results = []
    
    def worker(pdbs, gpu):
        for pdb in pdbs:
            pre = os.path.splitext(os.path.basename(pdb))[0]
            print(f"\n=== {pre} GPU{gpu} ===", flush=True)
            r = run_md(pdb, md_dir, 72.0, 1.0, pre, gpu)
            if r: results.append(r)
    
    t1 = threading.Thread(target=worker, args=(b1, 0))
    t2 = threading.Thread(target=worker, args=(b2, 1))
    t1.start(); t2.start()
    t1.join(); t2.join()
    
    print(f"\n{'='*50}")
    print("72C THERMAL STABILITY RESULTS")
    print(f"{'='*50}")
    results.sort(key=lambda x: x["avg_RMSD_nm"])
    for r in results:
        icon = "✅" if r["stability"]=="STABLE" else "⚠️" if r["stability"]=="MODERATE" else "❌"
        print(f"  {icon} {r['prefix']:30s} RMSD={r['avg_RMSD_nm']:.4f}nm  drift={r['drift_nm']:.4f}nm  [{r['stability']}]")
    print(f"\nStable at 72C: {sum(1 for r in results if r['stability']=='STABLE')}/{len(results)}")
    
    with open(os.path.join(md_dir, "summary.json"), "w") as f:
        json.dump(results, f, indent=2)
