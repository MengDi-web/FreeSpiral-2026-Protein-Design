# FreeSpiral - 2026 Synthetic Biology Innovation Competition

## Team FreeSpiral: AI-Guided GFP Design for Optimal Brightness and Thermal Stability

### Competition Overview

This repository contains our submission for the **2026 Protein Design Challenge**, part of the Synthetic Biology Innovation Competition. The goal is to design 6 GFP (Green Fluorescent Protein) variants (220-250 aa) that maximize:

**Score = Relative Brightness × Thermal Stability Retention**

where:
- **Relative Brightness** = F_initial / F_initial_WT (normalized to sfGFP)
- **Thermal Stability** = F_final / F_initial (after 72°C for 10 min)
- Only sequences with F_initial > 0.3 × F_initial_WT are eligible

### Repository Structure

```
├── README.md                 # This file
├── submission.csv            # Final 6 sequences for competition
├── docs/
│   ├── design_report.pdf     # Full design methodology document
│   └── analysis/
│       ├── brightness_analysis.ipynb  # Data exploration notebook
│       └── mutation_impact.ipynb      # Mutation effect analysis
├── scripts/
│   ├── design_sequences.py   # Sequence generation & validation
│   ├── validate.py           # Competition rule validation
│   └── exclusion_check.py    # Exclusion list validation
├── data/
│   ├── GFP_data.xlsx         # Provided brightness dataset
│   ├── Exclusion_List.csv    # Prohibited sequences
│   └── AAseqs of 5 GFP proteins.txt  # Reference sequences
└── results/
    └── designed_sequences.txt # Annotated designed sequences
```

### 6 Designed Sequences

| ID | Strategy | Backbone | Key Mutations |
|----|----------|----------|---------------|
| Seq1 | **Champion** (Balanced) | avGFP | S65T, S72A, Q80R, Y143G, S147P, Q157G, V163A, I167T, I171V, S202D |
| Seq2 | **BrightStar** (Max Brightness) | avGFP | S30R, S65T, S72A, Q80R, Y143G, S147P, Q157G, V163A, I167T, I171V |
| Seq3 | **BrightPlus** (sfGFP-based) | sfGFP | S72A, Y143G, S147P, Q157G, I167T, I171V |
| Seq4 | **Stable** (Max Stability) | sfGFP | L18E, S147P, Q157G, I167T, I171V |
| Seq5 | **Conservative** (Winner-proven) | avGFP | S30R, L64F, S65T, S72A, Q80R, S147P, V163A, I167T, I171V, S202D |
| Seq6 | **Evolved** (Diverse) | avGFP | S65T, S72A, Q80R, Y143G, Q157G, I167T, I171V, T203I |

### Design Methodology

Our approach integrates three complementary strategies:

#### 1. Data-Driven Rational Design
- Analyzed 51,715 avGFP mutation-brightness measurements from GFP_data.xlsx
- Identified single mutations with highest fold-improvement over WT (K157G: +2.48×, S174T: +1.94×, R72L: +1.88×)
- Cross-referenced with the 20 previous top-10 sequences (2024-2025) to validate mutation efficacy

#### 2. Literature-Backed Stability Engineering
- Integrated superfolder GFP mutations (S30R, Y39N, S65T, F99S, N105T, Y145F, M153T, V163A) as baseline
- Applied TGP-inspired surface charge optimization for thermal stability
- Incorporated loop-stabilizing proline substitutions (S147P) and core packing improvements (I167T, I171V)

#### 3. Knowledge-Based Design via LLM Agent
- Used Claude/GPT-5 (Codex) to execute an automated design pipeline
- Analyzed sequence-function relationships, literature review, and mutation compatibility
- Applied multi-strategy selection with cross-validation against exclusion list (135,414 sequences)

### Validation

All 6 sequences pass competition requirements:
- ✅ Length: 238-239 aa (within 220-250)
- ✅ Start with Methionine (M)
- ✅ Only standard 20 amino acids
- ✅ No stop codons
- ✅ Not in Exclusion_List.csv
- ✅ Diverse backbones (avGFP & sfGFP) for risk mitigation

### Environment Setup

```bash
# Python 3.9+ required
pip install openpyxl pandas pypdf
```

### Quick Start

```python
# Validate sequences against competition rules
python scripts/validate.py

# Check exclusion list
python scripts/exclusion_check.py

# Generate new designs
python scripts/design_sequences.py
```

### Acknowledgments

- Competition organizers for providing GFP_data.xlsx and reference materials
- Reference papers: Superfolder GFP (Pédelacq et al., 2006), TGP (Close et al., 2015)

### License

MIT
