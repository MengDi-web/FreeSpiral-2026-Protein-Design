#!/usr/bin/env python3
"""
Lightweight MD Simulation using OpenMM (no GROMACS needed)
模拟72°C热稳定性
"""
import os, sys, json
import numpy as np

def run_md_simulation(pdb_file, output_dir, temp_celsius=72.0, 
                      sim_time_ns=5.0, prefix="protein"):
    """
    Run MD simulation using OpenMM
    """
    try:
        from openmm import app, unit, Platform
        from openmm.app import PDBFile, Modeller, ForceField, Simulation
        from openmm import (
            LangevinIntegrator, CustomExternalForce, 
            MonteCarloBarostat, System
        )
    except ImportError:
        print("OpenMM not installed. Install with: pip install openmm pdbfixer")
        return None
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Temperature
    temp = temp_celsius + 273.15  # Kelvin
    
    try:
        # Load PDB
        pdb = app.PDBFile(pdb_file)
        
        # Force field
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        # Modeller
        modeller = Modeller(pdb.topology, pdb.positions)
        
        # Add solvent
        modeller.addSolvent(forcefield, padding=1.2*unit.nanometer)
        modeller.addHydrogens(forcefield)
        
        # System
        system = forcefield.createSystem(
            modeller.topology, 
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds
        )
        
        # Barostat for NPT
        system.addForce(MonteCarloBarostat(1.0*unit.atmosphere, temp, 25))
        
        # Integrator
        integrator = LangevinIntegrator(
            temp*unit.kelvin,
            1.0/unit.picosecond,
            0.002*unit.picoseconds
        )
        
        # Platform
        try:
            platform = Platform.getPlatformByName('CUDA')
            print("  Using CUDA platform")
        except:
            try:
                platform = Platform.getPlatformByName('OpenCL')
                print("  Using OpenCL platform")
            except:
                platform = Platform.getPlatformByName('CPU')
                print("  Using CPU platform")
        
        # Simulation
        simulation = Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        
        # Energy minimization
        print("  Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=5000)
        
        # Save minimized structure
        state = simulation.context.getState(getPositions=True)
        with open(os.path.join(output_dir, f"{prefix}_minimized.pdb"), 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Equilibrate at 30°C
        print("  Equilibrating at 30°C...")
        integrator.setTemperature(303.15*unit.kelvin)
        simulation.step(25000)  # 50ps
        
        # Heat to 72°C
        print(f"  Heating to {temp_celsius}°C...")
        n_steps_heat = 50000
        for i in range(10):
            t = 303.15 + (temp - 303.15) * (i+1) / 10
            integrator.setTemperature(t*unit.kelvin)
            simulation.step(n_steps_heat // 10)
        
        # Production at 72°C
        total_steps = int(sim_time_ns * 500)  # 500 steps/ns at 2fs
        n_save = max(1, total_steps // 100)  # Save 100 frames
        
        print(f"  Running production at {temp_celsius}°C ({sim_time_ns}ns)...")
        
        rmsd_data = []
        rmsf_data = {}
        
        # Reference positions for RMSD
        ref_state = simulation.context.getState(getPositions=True)
        ref_positions = ref_state.getPositions()
        
        for step in range(total_steps):
            simulation.step(1)
            
            if step % n_save == 0 or step == total_steps - 1:
                state = simulation.context.getState(getPositions=True, getEnergy=True)
                current_pos = state.getPositions()
                
                # Calculate RMSD
                rmsd = np.sqrt(np.mean([
                    (ref_positions[i][j].value_in_unit(unit.nanometer) - 
                     current_pos[i][j].value_in_unit(unit.nanometer))**2
                    for i in range(min(len(ref_positions), len(current_pos)))
                    for j in range(3)
                ]))
                rmsd_data.append(rmsd)
                
                # Track per-residue RMSF
                for i in range(min(len(ref_positions), len(current_pos))):
                    if i not in rmsf_data:
                        rmsf_data[i] = []
                    rmsf_data[i].append(np.sqrt(np.mean([
                        (ref_positions[i][j].value_in_unit(unit.nanometer) - 
                         current_pos[i][j].value_in_unit(unit.nanometer))**2
                        for j in range(3)
                    ])))
                
                if step % (n_save * 10) == 0:
                    ns_done = step * 0.002 / 1000  # ps -> ns
                    print(f"    {ns_done:.1f}ns / {sim_time_ns}ns (RMSD={rmsd:.3f}nm)")
        
        # Final structure
        state = simulation.context.getState(getPositions=True)
        with open(os.path.join(output_dir, f"{prefix}_final.pdb"), 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        
        # Results
        avg_rmsd = np.mean(rmsd_data[-min(50, len(rmsd_data)//2):])
        max_rmsd = np.max(rmsd_data)
        rmsd_drift = rmsd_data[-1] - rmsd_data[0] if len(rmsd_data) > 1 else 0
        
        result = {
            "prefix": prefix,
            "temperature_C": temp_celsius,
            "sim_time_ns": sim_time_ns,
            "avg_RMSD_nm": round(float(avg_rmsd), 3),
            "max_RMSD_nm": round(float(max_rmsd), 3),
            "RMSD_drift_nm": round(float(rmsd_drift), 3),
            "stability": "STABLE" if avg_rmsd < 0.3 and rmsd_drift < 0.1 else \
                        "MODERATE" if avg_rmsd < 0.5 else "UNSTABLE"
        }
        
        # Save results
        with open(os.path.join(output_dir, f"{prefix}_md_results.json"), 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"\n  ✅ {prefix}: RMSD={avg_rmsd:.3f}nm, max={max_rmsd:.3f}nm [{result['stability']}]")
        return result
    
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None

if __name__ == "__main__":
    import glob
    
    pdb_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "md_results"
    
    pdbs = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    
    if not pdbs:
        print(f"No PDB files found in {pdb_dir}")
        print("Usage: python run_openmm_md.py <pdb_dir> <output_dir>")
        sys.exit(1)
    
    print(f"Found {len(pdbs)} PDB files")
    
    all_results = []
    for pdb_file in pdbs:
        prefix = os.path.splitext(os.path.basename(pdb_file))[0]
        print(f"\n{'='*50}")
        print(f"Running MD: {prefix}")
        print(f"{'='*50}")
        
        result = run_md_simulation(
            pdb_file=pdb_file,
            output_dir=output_dir,
            temp_celsius=72.0,
            sim_time_ns=5.0,
            prefix=prefix
        )
        if result:
            all_results.append(result)
    
    # Summary
    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    for r in all_results:
        print(f"  {r['prefix']}: RMSD={r['avg_RMSD_nm']:.3f}nm [{r['stability']}]")
    
    stable = sum(1 for r in all_results if r['stability'] == 'STABLE')
    print(f"\nStable at 72°C: {stable}/{len(all_results)}")
    
    with open(os.path.join(output_dir, "md_summary.json"), 'w') as f:
        json.dump({"results": all_results, "stable_count": stable}, f, indent=2)
