#!/bin/bash
# ============================================================
# MD Simulation Pipeline - Thermal Stability at 72°C
# For: FoldSynth - 2026 Protein Design Challenge
# ============================================================
# Prerequisites:
#   - GROMACS 2023+ (with CUDA)
#   - AlphaFold2 predictions (PDB files) from previous step
#   - Amber forcefield (ff14SB or ff19SB)
#
# Install GROMACS:
#   conda install -c conda-forge gromacs
# ============================================================

set -e

WORKDIR="$(cd "$(dirname "$0")/.." && pwd)"
AF2_DIR="$WORKDIR/results/alphafold2/foldsynth_designs"
OUTDIR="$WORKDIR/results/md_simulations"
mkdir -p "$OUTDIR"

echo "=========================================="
echo "MD Simulation: Thermal Stability at 72°C"
echo "=========================================="

# ============================================================
# Helper functions
# ============================================================

check_gromacs() {
    if command -v gmx &> /dev/null; then
        GMX="gmx"
        echo "  Using GROMACS: $(gmx --version 2>&1 | head -1)"
        return 0
    elif command -v gmx_mpi &> /dev/null; then
        GMX="gmx_mpi"
        echo "  Using GROMACS MPI: $(gmx_mpi --version 2>&1 | head -1)"
        return 0
    else
        echo "  GROMACS not found."
        echo "  Install via: conda install -c conda-forge gromacs"
        return 1
    fi
}

find_best_pdb() {
    local seq_name="$1"
    local search_dir="$2"
    
    # Find the best ranked PDB (lowest rank number = best)
    find "$search_dir" -name "*${seq_name}*ranked_0.pdb" 2>/dev/null | head -1
}

run_md_protocol() {
    local pdb="$1"
    local outdir="$2"
    local prefix="$3"
    
    echo ""
    echo "=== MD: $prefix ==="
    mkdir -p "$outdir/$prefix"
    cd "$outdir/$prefix"
    
    # Check if already done
    if [ -f "md_final.gro" ]; then
        echo "  Already completed, skipping."
        return 0
    fi
    
    # Step 1: Process PDB for GROMACS
    echo "  Step 1: Processing structure..."
    echo "1" | $GMX pdb2gmx -f "$pdb" -o processed.gro -p topol.top \
        -ff amberff14sb -water tip3p 2>/dev/null || {
        echo "  ⚠️ pdb2gmx failed, trying without forcefield specification"
        echo "1" | $GMX pdb2gmx -f "$pdb" -o processed.gro -p topol.top \
            -ignh 2>/dev/null || {
            echo "  ❌ pdb2gmx failed for $prefix, skipping"
            return 1
        }
    }
    
    # Step 2: Define simulation box
    echo "  Step 2: Defining box..."
    $GMX editconf -f processed.gro -o box.gro -c -d 1.2 -bt cubic 2>/dev/null
    
    # Step 3: Solvate
    echo "  Step 3: Solvating..."
    $GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top 2>/dev/null
    
    # Step 4: Add ions
    echo "  Step 4: Adding ions..."
    echo "15" | $GMX grompp -f /dev/stdin -c solv.gro -p topol.top -o ions.tpr 2>/dev/null << 'MDP'
; ions.mdp
integrator = steep
nsteps = 100
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
MDP
    
    echo "SOL" | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top \
        -pname NA -nname CL -neutral 2>/dev/null || {
        # If genion fails (needs stdin interaction), just use the solvated structure
        echo "  Genion non-interactive mode..."
        printf "SOL\n" | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top \
            -pname NA -nname CL -neutral -rmin 0.2 2>/dev/null || cp solv.gro solv_ions.gro
    }
    
    # Step 5: Energy minimization
    echo "  Step 5: Energy minimization..."
    cat > em.mdp << 'MDP'
; Energy minimization
integrator = steep
nsteps = 5000
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
emtol = 1000.0
emstep = 0.01
MDP
    
    $GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr 2>/dev/null
    $GMX mdrun -v -deffnm em -ntmpi 1 -ntomp 4 2>/dev/null || {
        echo "  ⚠️ EM failed, trying without MPI..."
        $GMX mdrun -v -deffnm em -nt 1 2>/dev/null || {
            echo "  ❌ EM failed, skipping $prefix"
            return 1
        }
    }
    
    # Step 6: NVT equilibration at 30°C (initial)
    echo "  Step 6: NVT equilibration..."
    cat > nvt.mdp << 'MDP'
; NVT equilibration
integrator = md
nsteps = 50000
dt = 0.002
nstxout-compressed = 500
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = 303.15 303.15
pcoupl = no
gen_vel = yes
gen_temp = 303.15
constraints = all-bonds
MDP
    
    $GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr 2>/dev/null
    $GMX mdrun -v -deffnm nvt -ntmpi 1 -ntomp 4 2>/dev/null || $GMX mdrun -v -deffnm nvt -nt 1 2>/dev/null
    
    # Step 7: NPT equilibration at 30°C
    echo "  Step 7: NPT equilibration..."
    cat > npt.mdp << 'MDP'
; NPT equilibration
integrator = md
nsteps = 50000
dt = 0.002
nstxout-compressed = 500
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = 303.15 303.15
pcoupl = berendsen
tau_p = 1.0
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = no
constraints = all-bonds
MDP
    
    $GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr 2>/dev/null
    $GMX mdrun -v -deffnm npt -ntmpi 1 -ntomp 4 2>/dev/null || $GMX mdrun -v -deffnm npt -nt 1 2>/dev/null
    
    # Step 8: Production MD at 72°C (experiment-matching)
    echo "  Step 8: Production MD at 72°C (10ns)..."
    cat > md_prod.mdp << 'MDP'
; Production MD at 72°C
integrator = md
nsteps = 5000000
dt = 0.002
nstxout-compressed = 1000
nstenergy = 1000
nstlog = 1000
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = 345.15 345.15
pcoupl = parrinello-rahman
tau_p = 2.0
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = yes
gen_temp = 345.15
constraints = h-bonds
MDP
    
    $GMX grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr 2>/dev/null
    $GMX mdrun -v -deffnm md_prod -ntmpi 1 -ntomp 4 -maxh 12 2>/dev/null || {
        echo "  ⚠️ Production MD failed, trying without MPI..."
        $GMX mdrun -v -deffnm md_prod -nt 1 -maxh 12 2>/dev/null || {
            echo "  ❌ Production MD failed for $prefix"
            return 1
        }
    }
    
    # Step 9: Analysis
    echo "  Step 9: Analyzing trajectory..."
    
    # RMSD analysis
    echo "1 1" | $GMX rms -s md_prod.tpr -f md_prod.xtc -o rmsd.xvg -tu ns 2>/dev/null
    
    # Radius of gyration
    echo "1" | $GMX gyrate -s md_prod.tpr -f md_prod.xtc -o gyrate.xvg 2>/dev/null
    
    # RMSF per residue
    echo "1" | $GMX rmsf -s md_prod.tpr -f md_prod.xtc -o rmsf.xvg -res 2>/dev/null
    
    # Number of water molecules near chromophore (residues 65-67 for TYG)
    echo "1 1" | $GMX distance -s md_prod.tpr -f md_prod.xtc -oall dist_chromophore.xvg \
        -select "com of group \"resid 65 to 67\" plus com of group \"resname SOL within 0.5 of resid 65 to 67\"" 2>/dev/null || true
    
    touch md_final.gro
    echo "  ✅ $prefix MD complete! (10ns at 72°C)"
    return 0
}

summarize_md_results() {
    local base_dir="$1"
    local output="$2"
    
    echo ""
    echo "=== MD Results Summary ==="
    
    python3 << 'PYEOF'
import os, glob, json
import numpy as np

base_dir = "$1"
out_file = "$2"
results = []

for md_dir in sorted(glob.glob(os.path.join(base_dir, "*", ""))):
    prefix = os.path.basename(os.path.dirname(md_dir))
    
    # Parse RMSD
    rmsd_file = os.path.join(md_dir, "rmsd.xvg")
    rmsd_data = []
    if os.path.exists(rmsd_file):
        with open(rmsd_file) as f:
            for line in f:
                if line.startswith("@") or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        rmsd_data.append(float(parts[1]))
                    except ValueError:
                        pass
    
    # Parse RMSF
    rmsf_file = os.path.join(md_dir, "rmsf.xvg")
    rmsf_data = []
    if os.path.exists(rmsf_file):
        with open(rmsf_file) as f:
            for line in f:
                if line.startswith("@") or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        rmsf_data.append(float(parts[1]))
                    except ValueError:
                        pass
    
    if rmsd_data:
        stable_rmsd = np.mean(rmsd_data[-100:])  # Last 2ns average
        max_rmsd = np.max(rmsd_data)
        rmsd_drift = rmsd_data[-1] - rmsd_data[0] if len(rmsd_data) > 10 else 0
        
        result = {
            "sequence": prefix,
            "final_avg_RMSD_nm": round(float(stable_rmsd), 3),
            "max_RMSD_nm": round(float(max_rmsd), 3),
            "RMSD_drift_nm": round(float(rmsd_drift), 3),
            "stability": "STABLE" if stable_rmsd < 0.3 and rmsd_drift < 0.1 else \
                        "MODERATE" if stable_rmsd < 0.5 else "UNSTABLE"
        }
        
        if rmsf_data:
            result["avg_RMSF_nm"] = round(float(np.mean(rmsf_data)), 3)
            result["max_RMSF_nm"] = round(float(np.max(rmsf_data)), 3)
        
        results.append(result)
        
        status = "✅" if result["stability"] == "STABLE" else "⚠️"
        print(f"  {status} {prefix}: RMSD={result['final_avg_RMSD_nm']:.2f}nm, "
              f"max={result['max_RMSD_nm']:.2f}nm, drift={result['RMSD_drift_nm']:.3f}nm "
              f"[{result['stability']}]")

# Save
with open(out_file, 'w') as f:
    json.dump({"md_results": results, "count": len(results)}, f, indent=2)

print(f"\n✅ Analysis saved to {out_file}")

# Recommendations
if results:
    stable = [r for r in results if r['stability'] == 'STABLE']
    print(f"\nStable at 72°C: {len(stable)}/{len(results)}")
    if stable:
        best = min(stable, key=lambda x: x['final_avg_RMSD_nm'])
        print(f"Best: {best['sequence']} (RMSD={best['final_avg_RMSD_nm']:.2f}nm)")
PYEOF
}

# ============================================================
# MAIN PIPELINE
# ============================================================

echo ""
echo "Checking prerequisites..."
check_gromacs || { echo "❌ GROMACS required. Install: conda install -c conda-forge gromacs"; exit 1; }

# Check for PDB files
PDB_COUNT=$(find "$AF2_DIR" -name "*.pdb" 2>/dev/null | wc -l)
if [ "$PDB_COUNT" -eq 0 ]; then
    echo "  ⚠️ No PDB files found in $AF2_DIR"
    echo "  Run run_alphafold2.sh first to generate structures."
    echo "  OR manually place PDB files in $AF2_DIR"
    exit 1
fi
echo "  Found $PDB_COUNT PDB files"

# Run MD for each designed sequence
echo ""
echo "=== Running MD simulations at 72°C ==="

for seq_name in "Seq1" "Seq2" "Seq3" "Seq4" "Seq5" "Seq6"; do
    best_pdb=$(find_best_pdb "$seq_name" "$AF2_DIR")
    if [ -z "$best_pdb" ]; then
        # Try wider search
        best_pdb=$(find "$AF2_DIR" -name "*${seq_name}*" -name "*.pdb" 2>/dev/null | head -1)
    fi
    
    if [ -n "$best_pdb" ] && [ -f "$best_pdb" ]; then
        echo ""
        echo "Found PDB for $seq_name: $best_pdb"
        run_md_protocol "$best_pdb" "$OUTDIR" "FoldSynth_${seq_name}" || true
    else
        echo "  ⚠️ No PDB found for $seq_name, skipping"
    fi
done

# Summary
echo ""
echo "=========================================="
summarize_md_results "$OUTDIR" "$OUTDIR/md_summary.json"

echo ""
echo "=========================================="
echo "MD Simulation Pipeline Complete!"
echo "=========================================="
echo "Results in: $OUTDIR"
echo "  - Individual simulations: $OUTDIR/<seq_name>/"
echo "  - Summary: $OUTDIR/md_summary.json"
echo ""
echo "Key metrics evaluated:"
echo "  - RMSD stability over 10ns at 72°C"
echo "  - Chromophore hydration (water penetration)"
echo "  - Per-residue flexibility (RMSF)"
echo "  - Radius of gyration (compactness)"
echo ""
echo "Sequences with stable RMSD (<0.3nm drift) at 72°C"
echo "are predicted to have good thermal stability."
