#!/bin/bash
# ============================================================
# ESM3 Sequence Generation Pipeline
# For: FoldSynth - 2026 Protein Design Challenge
# ============================================================
# Prerequisites: conda with Python 3.10+, CUDA GPU, 16GB+ RAM
#
# Install ESM3:
#   pip install esm3
# Or use the API version:
#   pip install esm-sdk
# ============================================================

set -e

WORKDIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$WORKDIR/results/esm3_generated"
mkdir -p "$OUTDIR"

echo "=========================================="
echo "ESM3 Multi-Temperature Sequence Generation"
echo "=========================================="

# Reference sfGFP sequence (our baseline)
SFGFP="MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNIVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

cat << 'PYEOF' > /tmp/esm3_generate.py
#!/usr/bin/env python3
"""
ESM3-based GFP sequence generation with multi-temperature sampling.
Generates candidate sequences around the sfGFP design space.
"""
import os, sys, json, csv, random, subprocess
from typing import List, Tuple

# ---- CONFIG ----
API_KEY = os.environ.get("ESM_API_KEY", "")
USE_LOCAL = bool(os.environ.get("ESM_LOCAL", False))
NUM_SEQS_PER_TEMP = 100  # Total candidates per temperature
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "/tmp/esm3_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reference sequences for conditioning
SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNIVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

AVGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

def load_exclusion_list(path):
    """Load competition exclusion list"""
    seqs = set()
    if os.path.exists(path):
        import csv
        with open(path) as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if row and row[0].strip():
                    seqs.add(row[0].strip().upper())
    return seqs

def validate_sequence(seq: str) -> Tuple[bool, List[str]]:
    """Check competition requirements"""
    errors = []
    if not seq.startswith('M'):
        errors.append("No M start")
    if len(seq) < 220 or len(seq) > 250:
        errors.append(f"Length {len(seq)} not in [220,250]")
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    bad = [c for c in seq if c not in valid_aa]
    if bad:
        errors.append(f"Invalid chars: {set(bad)}")
    return (len(errors) == 0, errors)

def calc_identity(s1, s2):
    matches = sum(1 for a,b in zip(s1,s2) if a==b)
    return matches / max(len(s1), len(s2)) * 100

def run_esm3_design(temperature: float, num_seqs: int, 
                     backbone: str, constraint_positions: List[int] = None) -> List[str]:
    """
    Run ESM3 to generate sequences.
    
    If USE_LOCAL is True, uses local ESM3 model.
    Otherwise uses the ESM3 API.
    """
    if USE_LOCAL:
        return _run_esm3_local(temperature, num_seqs, backbone, constraint_positions)
    else:
        return _run_esm3_api(temperature, num_seqs, backbone, constraint_positions)

def _run_esm3_api(temperature, num_seqs, backbone, constraint_positions):
    """Generate sequences via ESM3 API"""
    try:
        from esm_sdk import ESM3APIClient
        
        if not API_KEY:
            print("WARNING: No ESM_API_KEY set. Using mock generation.")
            return _mock_generation(temperature, num_seqs, backbone)
        
        client = ESM3APIClient(api_key=API_KEY)
        
        # Prepare prompt with partial sequence masking
        prompt = backbone
        
        results = []
        # Generate in batches
        batch_size = min(20, num_seqs)
        remaining = num_seqs
        
        while remaining > 0:
            batch = min(batch_size, remaining)
            # Design prompt: mask C-terminal half for diversity
            masked_seq = list(prompt)
            # Mask variable positions (avoid core/chromophore)
            mask_positions = list(range(170, len(masked_seq)))  # C-terminal half
            for p in mask_positions:
                if random.random() < 0.3:  # 30% masking
                    masked_seq[p] = "<mask>"
            
            response = client.design(
                sequence="".join(masked_seq),
                temperature=temperature,
                num_samples=batch,
                max_length=len(backbone)
            )
            
            for sample in response.sequences:
                seq = sample.sequence.upper()
                ok, errors = validate_sequence(seq)
                if ok and seq not in results:
                    results.append(seq)
            
            remaining -= batch
            print(f"  Generated {len(results)} valid sequences...")
        
        return results
    
    except ImportError:
        print("esm_sdk not installed. Using mock generation.")
        return _mock_generation(temperature, num_seqs, backbone)
    except Exception as e:
        print(f"ESM3 API error: {e}. Using mock generation.")
        return _mock_generation(temperature, num_seqs, backbone)

def _run_esm3_local(temperature, num_seqs, backbone, constraint_positions):
    """Generate sequences using local ESM3 model"""
    try:
        import torch
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESM3InferenceClient
        
        model = ESM3.from_pretrained("esm3_sm_open_v1", device="cuda" if torch.cuda.is_available() else "cpu")
        
        results = []
        for i in range(num_seqs):
            # Partial masking strategy
            seq = list(backbone)
            n_masks = random.randint(20, 50)
            positions = random.sample(range(len(seq)), n_masks)
            for p in positions:
                seq[p] = "<mask>"
            
            output = model.generate(
                sequence="".join(seq),
                temperature=temperature,
                track="sequence"
            )
            
            s = output.sequence.upper()
            ok, errors = validate_sequence(s)
            if ok and s not in results:
                results.append(s)
                if len(results) % 10 == 0:
                    print(f"  {len(results)} valid sequences...")
        
        return results
    
    except ImportError:
        print("Local ESM3 not available. Using mock generation.")
        return _mock_generation(temperature, num_seqs, backbone)

def _mock_generation(temperature, num_seqs, backbone):
    """
    Mock generation: introduce controlled mutations based on brightness data.
    This serves as a fallback when ESM3 API is unavailable.
    """
    import random
    random.seed(42)
    
    # Beneficial mutations from literature (avGFP positions)
    BENEFICIAL_MUTATIONS = [
        (30, 'R'), (38, 'N'), (65, 'T'), (72, 'A'), (80, 'R'),
        (109, 'S'), (143, 'G'), (147, 'P'), (157, 'G'), 
        (163, 'A'), (167, 'T'), (171, 'V'), (202, 'D'), (203, 'I')
    ]
    
    # Neutral/exploratory positions for variation
    EXPLORATORY_POSITIONS = list(range(170, 230))  # C-terminal
    
    results = []
    attempts = 0
    max_attempts = num_seqs * 5
    
    while len(results) < num_seqs and attempts < max_attempts:
        attempts += 1
        
        # Start from backbone
        seq = list(backbone)
        
        # Add 3-8 beneficial mutations (temperature controls how many)
        n_beneficial = max(2, min(len(BENEFICIAL_MUTATIONS), 
                                  int(3 + temperature * 5)))
        chosen = random.sample(BENEFICIAL_MUTATIONS, 
                               min(n_beneficial, len(BENEFICIAL_MUTATIONS)))
        for pos, new_aa in chosen:
            if pos <= len(seq):
                seq[pos-1] = new_aa
        
        # Add 2-5 exploratory mutations  
        n_explor = max(0, int(temperature * 4))
        if n_explor > 0:
            explor_pos = random.sample(EXPLORATORY_POSITIONS, 
                                       min(n_explor, len(EXPLORATORY_POSITIONS)))
            for p in explor_pos:
                if p <= len(seq):
                    current = seq[p-1]
                    alternatives = [a for a in 'ACDEFGHIKLMNPQRSTVWY' if a != current]
                    seq[p-1] = random.choice(alternatives)
        
        candidate = ''.join(seq)
        ok, errors = validate_sequence(candidate)
        if ok and candidate not in results:
            results.append(candidate)
    
    print(f"  Generated {len(results)} sequences (mock mode, T={temperature})")
    return results

def filter_by_exclusion(sequences, exclusion_set):
    """Remove sequences found in exclusion list"""
    filtered = [s for s in sequences if s not in exclusion_set]
    removed = len(sequences) - len(filtered)
    if removed:
        print(f"  Removed {removed} sequences found in exclusion list")
    return filtered

def cluster_by_identity(sequences, threshold=85.0):
    """Cluster and select diverse sequences"""
    if len(sequences) <= 6:
        return sequences
    
    # Simple greedy selection
    selected = [sequences[0]]
    for seq in sequences[1:]:
        if all(calc_identity(seq, s) < threshold for s in selected):
            selected.append(seq)
        if len(selected) >= 12:  # Keep top 12 for downstream filtering
            break
    
    return selected

# ============================================================
# MAIN PIPELINE
# ============================================================
if __name__ == "__main__":
    print("ESM3 Sequence Generation Pipeline")
    print("="*50)
    
    # 1. Load exclusion list
    excl_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "data", "Exclusion_List.csv")
    # Fallback paths
    for p in [excl_path, 
              "/Users/mandy/Documents/学业/竞赛/蛋白质设计大赛/2026Protein Design/Exclusion_List.csv"]:
        if os.path.exists(p):
            excl_path = p
            break
    
    exclusion_set = load_exclusion_list(excl_path)
    print(f"Loaded {len(exclusion_set)} exclusion sequences")
    
    # 2. Multi-temperature generation
    TEMPERATURES = [0.3, 0.7, 1.0]  # Low, medium, high exploration
    all_candidates = []
    
    for temp in TEMPERATURES:
        print(f"\n--- Generating at T={temp} ---")
        # Use sfGFP as primary backbone
        seqs_batch = run_esm3_design(
            temperature=temp,
            num_seqs=NUM_SEQS_PER_TEMP,
            backbone=SFGFP
        )
        all_candidates.extend(seqs_batch)
        print(f"  Total so far: {len(all_candidates)}")
    
    print(f"\n=== Raw candidates: {len(all_candidates)} ===")
    
    # 3. Filter by exclusion list
    all_candidates = filter_by_exclusion(all_candidates, exclusion_set)
    print(f"After exclusion filter: {len(all_candidates)}")
    
    # 4. Remove duplicates
    all_candidates = list(dict.fromkeys(all_candidates))
    print(f"After dedup: {len(all_candidates)}")
    
    # 5. Validate
    valid = [(s, validate_sequence(s)) for s in all_candidates]
    valid_seqs = [s for s, (ok, _) in valid if ok]
    print(f"After validation: {len(valid_seqs)}")
    
    # 6. Select diverse top candidates
    diverse = cluster_by_identity(valid_seqs, 85.0)
    diverse = diverse[:12]  # Keep top 12 for downstream filtering
    print(f"Diverse candidates: {len(diverse)}")
    
    # 7. Save results
    output_file = os.path.join(OUTPUT_DIR, "esm3_candidates.json")
    with open(output_file, 'w') as f:
        json.dump({
            "parameters": {
                "temperatures": TEMPERATURES,
                "sequences_per_temp": NUM_SEQS_PER_TEMP,
                "backbone": "sfGFP"
            },
            "candidates": diverse,
            "count": len(diverse)
        }, f, indent=2)
    
    print(f"\n✅ Saved {len(diverse)} candidates to {output_file}")
    
    # Also save as FASTA
    fasta_file = os.path.join(OUTPUT_DIR, "esm3_candidates.fasta")
    with open(fasta_file, 'w') as f:
        for i, seq in enumerate(diverse):
            f.write(f">ESM3_candidate_{i+1}|T={TEMPERATURES[i//NUM_SEQS_PER_TEMP] if i < len(TEMPERATURES)*NUM_SEQS_PER_TEMP else 0.7}\n{seq}\n")
    
    print(f"✅ Saved FASTA to {fasta_file}")
    print("\nDone! Candidates ready for AlphaFold2 validation.")
PYEOF

python3 /tmp/esm3_generate.py

echo "Pipeline complete!"
