#!/usr/bin/env python3
# FreeSpiral - Competition Validation Script
# Validates submission.csv against ALL 2026 Protein Design rules

import sys, os, csv, argparse
from datetime import datetime

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED_HEADER = ["Team_Name", "Seq_ID", "Sequence"]
MIN_LEN, MAX_LEN = 220, 250

def check_seq(seq, exclusion_set):
    errors = []
    if not seq.startswith("M"):
        errors.append("Must start with M")
    length = len(seq)
    if length < MIN_LEN or length > MAX_LEN:
        errors.append(f"Length {length} outside [{MIN_LEN}-{MAX_LEN}]")
    bad = sorted(set(c for c in seq if c not in STANDARD_AA))
    if bad:
        errors.append(f"Invalid chars: {bad}")
    if "*" in seq:
        errors.append("Stop codon found")
    if exclusion_set and seq in exclusion_set:
        errors.append("IN EXCLUSION LIST")
    return errors, length

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="submission.csv")
    p.add_argument("--exclusion", default=None)
    p.add_argument("--template", action="store_true")
    args = p.parse_args()

    if args.template:
        print("CSV: Team_Name,Seq_ID,Sequence")
        print("FreeSpiral,1,MSKGEELFTG...")
        return

    # Read CSV
    try:
        with open(args.csv) as f:
            reader = csv.DictReader(f)
            seqs = [row for row in reader if all(k in row and row[k].strip() for k in REQUIRED_HEADER)]
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"Read {len(seqs)} sequences")

    # Load exclusion list
    exclusion_set = None
    if args.exclusion:
        excl_path = args.exclusion
    else:
        for p in ["../2026Protein Design/Exclusion_List.csv",
                  "/Users/mandy/Documents/学业/竞赛/蛋白质设计大赛/2026Protein Design/Exclusion_List.csv",
                  "data/Exclusion_List.csv"]:
            if os.path.exists(p):
                excl_path = p
                break
        else:
            excl_path = None

    if excl_path and os.path.exists(excl_path):
        with open(excl_path) as f:
            reader = csv.reader(f)
            next(reader, None)
            exclusion_set = set(row[0].strip() for row in reader if row and row[0].strip())
        print(f"Exclusion list: {len(exclusion_set):,} sequences")
    else:
        print("Skipping exclusion check (file not found)")

    # Validate
    all_ok = True
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print()
    print("="*72)
    print("  FreeSpiral - COMPETITION SEQUENCE VALIDATION")
    print(f"  {now}")
    print("="*72)

    for row in seqs:
        sid = row["Seq_ID"].strip()
        seq = row["Sequence"].strip()
        errors, length = check_seq(seq, exclusion_set)
        status = "PASS" if not errors else "FAIL"
        print(f"  [{status}] Seq{sid} ({length} aa)")
        for e in errors:
            print(f"         {e}")
            all_ok = False

    print()
    print("="*72)
    print(f"  Overall: {'ALL PASS' if all_ok else 'FAILURES FOUND'}")
    print("="*72)
    sys.exit(0 if all_ok else 1)

if __name__ == "__main__":
    main()
