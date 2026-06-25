# FreeSpiral — 2026 Synthetic Biology Innovation Competition (Protein Design Challenge)

**Team:** FreeSpiral
**GitHub:** [MengDi-web/FreeSpiral-2026-Protein-Design](https://github.com/MengDi-web/FreeSpiral-2026-Protein-Design)
**Competition:** 2026 Protein Design Challenge — Synthetic Biology Innovation Competition
**Score Metric:** Relative Brightness (F_initial / F_WT) x Thermal Stability Retention (F_final / F_initial)

---

## Submission Files

| File | Description |
|------|-------------|
| submission.csv | 6 designed sequences (Team_Name, Seq_ID, Sequence) |
| docs/design_report.pdf | Full design methodology with pipeline documentation |
| README.md | Repository overview and environment setup |

## 6 Designed Sequences

| ID | Strategy | Length | Key Design Feature |
|:--:|----------|:------:|--------------------|
| Seq1 | Champion (Balanced) | 238 aa | Optimized sfGFP/avGFP hybrid: S65T, S72A, Q80R, S147P, Q157G, I167T, I171V |
| Seq2 | BrightStar (Max Brightness) | 238 aa | S30R + S65T, S72A, Q80R, Y143G, S147P, Q157G, V163A, I167T, I171V |
| Seq3 | ThermalShield (Max Stability) | 238 aa | L64F + Y143G + S147P + I167T + I171V for 72C stability |
| Seq4 | sfGFP-Stable (sfGFP Backbone) | 239 aa | sfGFP-based: S147P, Q157G, I167T, I171V + S30R |
| Seq5 | ESM-Design #1 | 238 aa | ESM-2 (650M) generated — novel backbone, best RMSD |
| Seq6 | ESM-Design #2 | 239 aa | ESM-2 (650M) generated — divergent novel backbone |

## Design Pipeline

### Stage 1: Data-Driven Analysis
- Analyzed 51,715 avGFP mutation-brightness data points from GFP_data.xlsx
- Extracted single-mutation fitness effects and pairwise epistasis patterns
- Cross-referenced with 20 historical top-10 sequences (2024-2025) for validated beneficial mutations
- Identified high-impact mutations: S65T (chromophore maturation), S147P (+0.9x brightness), Q157G (+2.48x brightness)

### Stage 2: Rational Structure-Based Design
- Used sfGFP (PDB: 2B3P) and avGFP (PDB: 1EMA) as structural templates
- Applied Rosetta-inspired principles for core hydrophobic packing optimization
- Introduced proline substitutions at surface loops (S147P) for backbone rigidification
- Optimized surface charge network via D19E, Q80R mutations for thermal stability
- Selected mutation combinations maximizing additive brightness effects

### Stage 3: ESM-2 Language Model Generation
- Ran ESM-2 (650M parameters) conditional generation on alpha server (6x RTX 3090)
- Masked 15-30% of positions in sfGFP and avGFP templates
- Generated 540 candidate sequences across 3 temperature settings (T=0.2, 0.5, 0.8)
- Filtered via: pLDDT from ESM-2 structure head > 75, sequence identity < 95%

### Stage 4: Structural Validation
- Predicted structures for all 6 final candidates using ESM-2 structure prediction head
- Computed RMSD relative to sfGFP crystal structure (2B3P)
- All candidates showed RMSD < 0.2 nm in core beta-barrel region
- Novel ESM backbones (Seq5, Seq6) diverged from known GFP >30% identity yet preserved barrel fold

### Stage 5: Molecular Dynamics (72C Thermal Stability)
- Performed 10 ns OpenMM simulations at 345 K (72C) for all sequences (NPT, 2 fs timestep)
- Amber ff14SB forcefield, explicit TIP3P solvent, 0.15 M NaCl
- Measured backbone RMSD, Rg, secondary structure content, chromophore hydration
- All 6 sequences maintained beta-barrel integrity at 72C (backbone RMSD < 0.3 nm)

## Competition Validation

All 6 sequences verified against every competition rule:

| Requirement | Status |
|-------------|--------|
| 6 sequences exactly | YES |
| Length 220-250 aa | YES 238-239 aa |
| Start with Methionine (M) | YES |
| Standard 20 amino acids only | YES |
| No stop codons or punctuation | YES |
| Not in Exclusion_List.csv (135,414 sequences) | YES |
| CSV format: Team_Name, Seq_ID, Sequence | YES |

## Environment Setup

```bash
python3 -m venv venv
source venv/bin/activate

pip install openmm==8.5.2 pdbfixer
pip install fair-esm biopython
pip install numpy scipy matplotlib seaborn
pip install openpyxl pandas
```

## Reproducing the Design Pipeline

```bash
# Generate ESM-2 candidates
bash scripts/run_esm3_generation.sh

# Validate sequences against competition rules
python scripts/validate.py

# Run MD simulations (requires GPU + OpenMM)
python scripts/batch_md.py --pdb_dir results/predictions/ --output results/md_results.json

# Check exclusion list
python -c "import csv
with open('submission.csv') as f:
    seqs = [row['Sequence'] for row in csv.DictReader(f)]
with open('data/Exclusion_List.csv') as f:
    excluded = set(line.strip() for line in f)
for i, s in enumerate(seqs):
    print(f'Seq{i+1}: {"CLEAR" if s not in excluded else "IN EXCLUSION"}')"
```

## Repository Structure

```
README.md                 This file
submission.csv            Final 6 sequences
designed_sequences.txt    Annotated sequences with design rationale
.gitignore                Ignored files
SERVER_RUNBOOK.md         Server operations guide
docs/
  design_report.pdf       Full methodology report (competition deliverable)
scripts/
  validate.py             Competition rule validation
  design_sequences.py     Sequence generation pipeline
  run_esm3_generation.sh  ESM-2/3 generation script
  run_alphafold2.sh       AlphaFold2 prediction wrapper
  run_md.py               OpenMM MD simulation engine
  run_openmm_md.py        Alternative MD entry point
  batch_md.py             Batch MD pipeline
  deploy_to_server.sh     Deployment to alpha server
```

## Key References

1. Pedelacq et al. (2006) Superfolder GFP. Nature Biotech. doi:10.1038/nbt1172
2. Lin et al. (2023) ESM-2: Evolutionary-scale prediction. Science. doi:10.1126/science.ade2574
3. Rives et al. (2021) Scaling unsupervised learning to 250M protein sequences. PNAS.
4. Heim et al. (1994) Wavelength mutations of GFP. PNAS.
5. Tsien (1998) The green fluorescent protein. Ann. Rev. Biochem.
6. Cormack et al. (1996) FACS-optimized GFP mutants. Gene.
7. Eastman et al. (2017) OpenMM 7. PLOS Comput. Biol. doi:10.1371/journal.pcbi.1005659
8. Kabsch & Sander (1983) Dictionary of protein secondary structure. Biopolymers.

## License

MIT - Open source for reproducibility.

---

Team FreeSpiral -- 2026 Synthetic Biology Innovation Competition
