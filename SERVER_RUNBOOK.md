# Server Runbook — FoldSynth Pipeline on Alpha Server

## Server Access

- **Host:** 172.20.73.14 (alpha)
- **User:** mengd
- **Home:** /mnt/home/mengd/
- **GPUs:** 6x NVIDIA RTX 3090 (24 GB each)
- **OS:** Ubuntu 22.04, CUDA 13.0

## Final Pipeline Executed

### 1. ESM-2 Sequence Generation
Generated 540 candidate sequences from sfGFP/avGFP templates using ESM-2 (650M).

### 2. Structure Prediction (ESM-2 Structure Head)
Predicted 3D structures for all candidates. Selected top 6 by pLDDT and RMSD.

### 3. MD Simulation (72C Thermal Stability)
Ran OpenMM 8.5.2 NPT simulations at 345 K (72C) for 10 ns per sequence.

## Important Notes

- Do NOT use GPUs 0-1 (reserved for other users). Pipeline uses GPUs 2-5.
- All project files in /mnt/home/mengd/FoldSynth/
- Virtual environment: /mnt/home/mengd/FoldSynth/.venv/
- Conda base available alongside venv

## Installed Tools

- Python 3.10
- OpenMM 8.5.2 + PDBFixer
- PyTorch 2.x + CUDA
- fair-esm (ESM-2)
- ColabFold
- NumPy, SciPy, Biopython

## Directory Structure on Server

```
/mnt/home/mengd/FoldSynth/
  submission.csv           Final sequences
  designed_sequences.fasta  FASTA format for ESM-2
  scripts/                 Pipeline scripts
  results/
    esm_candidates/        540 ESM-2 generated sequences
    predictions/           Predicted structures
    md_results/            MD simulation trajectories + logs
    batch_md_results.json  Batch MD stability results
```

## Quick Commands

```bash
# SSH
ssh mengd@172.20.73.14

# Activate environment
cd ~/FoldSynth
source .venv/bin/activate

# Check GPU availability
nvidia-smi

# Run MD pipeline
python scripts/batch_md.py --pdb_dir results/predictions/ --output results/md_results.json

# List all MD results
cat results/batch_md_results.json
```

## Competition Deliverables

Updated from server to local MacBook and pushed to GitHub.
