#!/bin/bash
# ============================================================
# AlphaFold2 Structure Prediction Pipeline
# For: FoldSynth - 2026 Protein Design Challenge
# ============================================================
# Prerequisites:
#   - AlphaFold2 installed (or ColabFold via conda)
#   - CUDA GPU recommended
#
# Install ColabFold (lighter alternative):
#   conda install -c bioconda -c conda-forge colabfold
# Or use docker for AlphaFold2
# ============================================================

set -e

WORKDIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$WORKDIR/results/alphafold2"
mkdir -p "$OUTDIR" "$WORKDIR/results/af2_predictions"
FASTA_DIR="$WORKDIR/results/esm3_generated"

echo "=========================================="
echo "AlphaFold2/ColabFold Structure Prediction"
echo "=========================================="

# Function to run ColabFold (if available)
run_colabfold() {
    local fasta="$1"
    local outdir="$2"
    local prefix="$3"
    
    echo "  Running ColabFold on $prefix..."
    
    colabfold_batch \
        --num-recycle 3 \
        --model-type alphafold2_multimer_v3 \
        --num-models 5 \
        --use-gpu-relax \
        --amber \
        --templates \
        --pdb-templates /path/to/pdb70 2>/dev/null \
        "$fasta" "$outdir" || {
        echo "  ColabFold not found, trying local AlphaFold2..."
        return 1
    }
    return 0
}

# Function to run AlphaFold2 directly
run_alphafold2() {
    local fasta="$1"  
    local outdir="$2"
    local prefix="$3"
    
    echo "  Running AlphaFold2 on $prefix..."
    
    # Check if alphafold is available
    if command -v alphafold &> /dev/null; then
        alphafold \
            --fasta_paths="$fasta" \
            --output_dir="$outdir" \
            --model_names=model_1,model_2,model_3,model_4,model_5 \
            --max_template_date=2026-01-01 \
            --db_preset=reduced_dbs \
            --use_gpu_relax=True
    elif command -v run_alphafold.sh &> /dev/null; then
        run_alphafold.sh \
            --fasta_paths="$fasta" \
            --output_dir="$outdir" \
            --model_names=model_1,model_2,model_3,model_4,model_5
    else
        echo "  Neither AlphaFold2 nor ColabFold found."
        echo "  Install via: conda install -c bioconda colabfold"
        echo "  Or download AlphaFold2 from https://github.com/google-deepmind/alphafold"
        return 1
    fi
}

# Function to analyze predictions
analyze_predictions() {
    local pdb_dir="$1"
    local output_file="$2"
    
    echo "  Analyzing predictions..."
    
    python3 << 'PYEOF' 
import os, glob, json
import numpy as np

pdb_dir = "$1" if len("$1") > 0 else "."
out_file = "$2" if len("$2") > 0 else "analysis.json"

results = []

# Parse pLDDT scores from PDB files (B-factor column)
for pdb_file in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
    base = os.path.basename(pdb_file)
    plddt_scores = []
    
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15] == "CA":
                try:
                    plddt = float(line[60:66].strip())
                    plddt_scores.append(plddt)
                except ValueError:
                    pass
    
    if plddt_scores:
        avg_plddt = np.mean(plddt_scores)
        min_plddt = np.min(plddt_scores)
        
        # Chromophore region (around residues 60-70 for TYG)
        chromo_region = slice(55, 75) if len(plddt_scores) > 75 else slice(0, len(plddt_scores))
        chromo_plddt = np.mean(plddt_scores[chromo_region]) if len(plddt_scores) > chromo_region.stop else avg_plddt
        
        results.append({
            "pdb": base,
            "avg_pLDDT": round(avg_plddt, 2),
            "chromophore_pLDDT": round(chromo_plddt, 2),
            "min_pLDDT": round(min_plddt, 2),
            "confidence": "HIGH" if avg_plddt > 90 else ("GOOD" if avg_plddt > 80 else "MODERATE")
        })
        
        status = "✅" if avg_plddt > 80 else "⚠️"
        print(f"  {status} {base}: avg_pLDDT={avg_plddt:.1f}, chromophore={chromo_plddt:.1f}")

# Save analysis
with open(out_file, 'w') as f:
    json.dump({"predictions": results, "count": len(results)}, f, indent=2)

# Summary
if results:
    best = max(results, key=lambda x: x['avg_pLDDT'])
    print(f"\n  Best: {best['pdb']} (avg_pLDDT={best['avg_pLDDT']})")
    high_conf = sum(1 for r in results if r['avg_pLDDT'] > 80)
    print(f"  High confidence (pLDDT>80): {high_conf}/{len(results)}")
    print(f"  ✅ Analysis saved to {out_file}")
PYEOF
}

# ============================================================
# MAIN PIPELINE
# ============================================================

# Step 1: Predict our FoldSynth designed sequences
echo ""
echo "Step 1: Predicting designed sequences..."
DESIGNED_FASTA="/tmp/designed_seqs.fasta"

python3 << 'PYEOF'
import csv, os

# Read the submission CSV and create FASTA
with open('$WORKDIR/submission.csv') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

with open('$DESIGNED_FASTA', 'w') as f:
    for row in rows:
        f.write(f">{row['Team_Name']}_{row['Seq_ID']}\n{row['Sequence']}\n")

print(f"Created FASTA with {len(rows)} sequences")
PYEOF

mkdir -p "$OUTDIR/foldsynth_designs"
if command -v colabfold_batch &> /dev/null; then
    run_colabfold "$DESIGNED_FASTA" "$OUTDIR/foldsynth_designs" "foldsynth"
elif command -v alphafold &> /dev/null || command -v run_alphafold.sh &> /dev/null; then
    run_alphafold2 "$DESIGNED_FASTA" "$OUTDIR/foldsynth_designs" "foldsynth"
else
    echo "⚠️ AlphaFold2 not found. Installing ColabFold..."
    pip install colabfold 2>/dev/null || conda install -y -c bioconda colabfold 2>/dev/null
    if command -v colabfold_batch &> /dev/null; then
        run_colabfold "$DESIGNED_FASTA" "$OUTDIR/foldsynth_designs" "foldsynth"
    else
        echo "❌ Could not install AlphaFold2/ColabFold."
        echo "   Install manually: conda install -c bioconda colabfold"
        echo "   Then re-run this script."
        exit 1
    fi
fi

# Step 2: Predict ESM3-generated candidates (if available)
if [ -f "$FASTA_DIR/esm3_candidates.fasta" ]; then
    echo ""
    echo "Step 2: Predicting ESM3 candidates..."
    mkdir -p "$OUTDIR/esm3_candidates"
    
    # Split into batches of 4 for GPU memory
    python3 << 'PYEOF'
import os
fasta = "$FASTA_DIR/esm3_candidates.fasta"
outdir = "$OUTDIR/esm3_candidates"

with open(fasta) as f:
    content = f.read()

seqs = content.strip().split('>')[1:]
batch_size = 4
batch = []
batch_num = 1

for seq_entry in seqs:
    lines = seq_entry.strip().split('\n')
    header = lines[0]
    sequence = ''.join(lines[1:])
    batch.append((header, sequence))
    
    if len(batch) >= batch_size:
        with open(f"{outdir}/batch_{batch_num}.fasta", 'w') as f:
            for h, s in batch:
                f.write(f">{h}\n{s}\n")
        batch = []
        batch_num += 1

# Remaining
if batch:
    with open(f"{outdir}/batch_{batch_num}.fasta", 'w') as f:
        for h, s in batch:
            f.write(f">{h}\n{s}\n")

print(f"Created {batch_num} batches for prediction")
PYEOF

    for batch_fa in "$OUTDIR/esm3_candidates"/batch_*.fasta; do
        bn=$(basename "$batch_fa" .fasta)
        mkdir -p "$OUTDIR/esm3_candidates/$bn"
        run_colabfold "$batch_fa" "$OUTDIR/esm3_candidates/$bn" "$bn" || \
        run_alphafold2 "$batch_fa" "$OUTDIR/esm3_candidates/$bn" "$bn" || \
            echo "  Skipping $bn (no AF2 available)"
    done
fi

# Step 3: Analyze all predictions
echo ""
echo "Step 3: Analyzing predictions..."
analyze_predictions "$OUTDIR/foldsynth_designs" "$OUTDIR/foldsynth_analysis.json"

if [ -d "$OUTDIR/esm3_candidates" ]; then
    analyze_predictions "$OUTDIR/esm3_candidates" "$OUTDIR/esm3_analysis.json"
fi

echo ""
echo "=========================================="
echo "AlphaFold2 Pipeline Complete!"
echo "=========================================="
echo "Results in: $OUTDIR"
echo "  - FoldSynth designs: $OUTDIR/foldsynth_designs/"
echo "  - ESM3 candidates: $OUTDIR/esm3_candidates/ (if generated)"
echo "  - Analysis: $OUTDIR/foldsynth_analysis.json"

# Print best scores
if [ -f "$OUTDIR/foldsynth_analysis.json" ]; then
    python3 -c "
import json
with open('$OUTDIR/foldsynth_analysis.json') as f:
    data = json.load(f)
print('\n=== FoldSynth Design pLDDT Scores ===')
for pred in sorted(data['predictions'], key=lambda x: x['avg_pLDDT'], reverse=True):
    print(f\"  {pred['pdb']}: avg={pred['avg_pLDDT']}, chromo={pred['chromophore_pLDDT']} [{pred['confidence']}]\")
"
fi
