#!/usr/bin/env python3
"""
Batch MD validation using thread_final.py's exact working pipeline.
Processes all sequences sequentially, avoids all threading issues.
"""
import os, sys, csv, json, math

WD = "/mnt/home/mengd/FoldSynth"
os.environ["CUDA_VISIBLE_DEVICES"] = "2"

AA_NAMES = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HID',
            'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
            'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
AA_ATOMS = {
    'ALA':['N','CA','C','O','CB'],'ARG':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'],
    'ASN':['N','CA','C','O','CB','CG','OD1','ND2'],'ASP':['N','CA','C','O','CB','CG','OD1','OD2'],
    'CYS':['N','CA','C','O','CB','SG'],'GLN':['N','CA','C','O','CB','CG','CD','OE1','NE2'],
    'GLU':['N','CA','C','O','CB','CG','CD','OE1','OE2'],'GLY':['N','CA','C','O'],
    'HID':['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2'],
    'ILE':['N','CA','C','O','CB','CG1','CG2','CD1'],'LEU':['N','CA','C','O','CB','CG','CD1','CD2'],
    'LYS':['N','CA','C','O','CB','CG','CD','CE','NZ'],'MET':['N','CA','C','O','CB','CG','SD','CE'],
    'PHE':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ'],
    'PRO':['N','CA','C','O','CB','CG','CD'],'SER':['N','CA','C','O','CB','OG'],
    'THR':['N','CA','C','O','CB','OG1','CG2'],
    'TRP':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],
    'TYR':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH'],
    'VAL':['N','CA','C','O','CB','CG1','CG2'],
}

# Parse 2B3P template
template = {}
with open(os.path.join(WD, "designs", "2B3P.pdb")) as f:
    for line in f:
        if line.startswith("ATOM") and line[21].strip() in ("A",""):
            rnum = int(line[22:26]); aname = line[12:16].strip()
            x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])
            if rnum not in template: template[rnum] = {}
            if aname not in template[rnum]: template[rnum][aname]=(x,y,z)
valid = sorted([r for r in template if "N" in template[r] and "CA" in template[r] and "C" in template[r]])

def run_md(seq, name):
    """Thread_final.py's exact pipeline for one sequence"""
    from pdbfixer import PDBFixer
    from openmm import app, unit
    from openmm.app import PDBFile, ForceField, Modeller, Simulation
    from openmm import LangevinIntegrator, MonteCarloBarostat, Platform
    import tempfile, math
    
    # Truncate to template length
    if len(seq) > len(valid):
        seq = seq[-len(valid):]
    seq_trunc = seq
    
    # Create PDB (EXACT format from thread_final.py)
    lines = ["CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1","MODEL        1"]
    an = 1
    for i, aa in enumerate(seq_trunc):
        rn = AA_NAMES.get(aa, "ALA")
        t = template[valid[i]]
        for nm in ["N","CA","C","O","CB"]:
            if nm in t and not (nm == "CB" and rn == "GLY"):
                x,y,z = t[nm]
                el = {"N":"N","O":"O"}.get(nm,"C")
                lines.append(f"ATOM  {an:5d} {nm:<4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
                an += 1
    lines += ["TER","ENDMDL"]
    pdb_content = "\n".join(lines)
    
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content); pdb_path = f.name
        
        fixer = PDBFixer(pdb_path)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        tmp_pdb = os.path.join(WD, "designs", "_tmp.pdb")
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(tmp_pdb, "w"))
        
        pdb2 = PDBFile(tmp_pdb)
        ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        mod = Modeller(pdb2.topology, pdb2.positions)
        mod.addHydrogens(ff)
        mod.addSolvent(ff, padding=0.5*unit.nanometer)
        
        system = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)
        system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, 300*unit.kelvin, 25))
        
        plat = Platform.getPlatformByName("OpenCL")
        integrator = LangevinIntegrator(303.15*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        sim = Simulation(mod.topology, system, integrator, plat, {"OpenCLPrecision":"mixed"})
        sim.context.setPositions(mod.positions)
        
        sim.minimizeEnergy(maxIterations=2000)
        sim.step(500)
        integrator.setTemperature(345.15*unit.kelvin)
        
        ref = sim.context.getState(getPositions=True).getPositions()
        rmsds = []
        for _ in range(10):
            sim.step(25000)
            cur = sim.context.getState(getPositions=True).getPositions()
            n = min(len(ref), len(cur))
            d = sum((ref[i][j].value_in_unit(unit.nanometer)-cur[i][j].value_in_unit(unit.nanometer))**2
                    for i in range(n) for j in range(3))
            rmsds.append(math.sqrt(d/(n*3)))
        
        os.unlink(pdb_path)
        final = rmsds[-1]
        print(f"  {name}: RMSD={final:.4f}nm", flush=True)
        del sim, integrator, system, mod, fixer
        import gc; gc.collect()
        return {"name":name, "rmsd":round(final,4), "status":"OK"}
    except Exception as e:
        print(f"  {name}: ERROR: {str(e)[:80]}", flush=True)
        return {"name":name, "rmsd":999, "status":str(e)[:80]}

if __name__ == "__main__":
    # Read our 6 designs
    seqs = []
    with open(os.path.join(WD, "submission.csv")) as f:
        for row in csv.DictReader(f):
            seqs.append((f"{row['Seq_ID']}_FreeSpiral", row['Sequence']))
    
    # Read top 10 ESM candidates
    esm_fasta = os.path.join(WD, "results", "esm_candidates", "esm_candidates.fasta")
    if os.path.exists(esm_fasta):
        with open(esm_fasta) as f:
            name, seq = "", ""
            for line in f:
                if line.startswith(">"):
                    if name and seq: seqs.append((name, seq))
                    name = line[1:30].strip(); seq = ""
                else: seq += line.strip()
            if name and seq: seqs.append((name, seq))
    
    print(f"Total sequences: {len(seqs)}")
    
    results = []
    for name, seq in seqs[:16]:  # 6 designs + 10 ESM
        r = run_md(seq, name)
        results.append(r)
        
        with open(os.path.join(WD, "results", "batch_md_results.json"), "w") as f:
            json.dump(sorted(results, key=lambda x: x["rmsd"]), f, indent=2)
    
    # Final ranking
    print(f"\n{'='*50}")
    print("BATCH MD RESULTS")
    print(f"{'='*50}")
    results.sort(key=lambda x: x["rmsd"])
    for r in results:
        tag = "✅" if r["status"]=="OK" else "❌"
        print(f"  {tag} {r['name'][:40]:40s} RMSD={r['rmsd']:.4f}")
