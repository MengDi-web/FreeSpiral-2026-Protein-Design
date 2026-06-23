#!/usr/bin/env python3
"""Validate designed sequences against competition rules"""
import csv, sys

SEQUENCES = []
with open('../submission.csv') as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        SEQUENCES.append((row[0], row[1], row[2]))

STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')

print("="*60)
print("COMPETITION RULE VALIDATION")
print("="*60)

all_pass = True
for team, seq_id, seq in SEQUENCES:
    errors = []
    
    # Rule 1: Must start with M
    if seq[0] != 'M':
        errors.append("❌ Does not start with Methionine (M)")
    
    # Rule 2: Length 220-250
    if len(seq) < 220 or len(seq) > 250:
        errors.append(f"❌ Length {len(seq)} not in [220, 250]")
    
    # Rule 3: Only standard amino acids
    invalid = [c for c in seq if c not in STANDARD_AA]
    if invalid:
        errors.append(f"❌ Invalid characters: {set(invalid)}")
    
    # Rule 4: No stop codons
    if '*' in seq:
        errors.append("❌ Contains stop codon (*)")
    
    print(f"\n{team} - {seq_id} (len={len(seq)})")
    if errors:
        for e in errors:
            print(f"  {e}")
        all_pass = False
    else:
        print(f"  ✅ All rules passed")

print(f"\n{'='*60}")
print(f"OVERALL: {'ALL PASS ✅' if all_pass else 'SOME FAILED ❌'}")
print(f"{'='*60}")

# Exclusion list check
print("\nChecking Exclusion_List.csv...")
try:
    with open('/Users/mandy/Documents/学业/竞赛/蛋白质设计大赛/2026Protein Design/Exclusion_List.csv') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        excl = {row[0].strip() for row in reader if row}
    
    print(f"  Exclusion list contains {len(excl)} sequences")
    for team, seq_id, seq in SEQUENCES:
        if seq in excl:
            print(f"  ❌ {seq_id}: FOUND IN EXCLUSION LIST")
            all_pass = False
        else:
            print(f"  ✅ {seq_id}: Not in exclusion list")
except FileNotFoundError:
    print("  ⚠️  Exclusion_List.csv not found, skipping")

sys.exit(0 if all_pass else 1)
