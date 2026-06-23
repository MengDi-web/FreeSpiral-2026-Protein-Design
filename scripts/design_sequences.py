#!/usr/bin/env python3
"""
FoldSynth Sequence Design Pipeline
AI-guided rational design of GFP variants for competition
"""
import csv

# Reference sequences
avGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
sfGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNIVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

def mutate(seq, changes):
    """Apply (pos, new_aa) mutations. pos is 1-indexed."""
    s = list(seq)
    for pos, new in changes:
        s[pos-1] = new
    return ''.join(s)

def validate(seq):
    """Check competition requirements"""
    return (seq.startswith('M') and 220<=len(seq)<=250 and 
            all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq))

# Design strategies with rationale
DESIGN_STRATEGIES = {
    "Seq1_avGFP_Champion": {
        "backbone": avGFP,
        "desc": "Balanced brightness+stability",
        "mutations": [
            (65,'T'), (72,'A'), (80,'R'), (143,'G'), (147,'P'),
            (157,'G'), (163,'A'), (167,'T'), (171,'V'), (202,'D')
        ],
        "rationale": "Top winner mutations + brightness enhancers"
    },
    "Seq2_avGFP_BrightStar": {
        "backbone": avGFP,
        "desc": "Maximum brightness",
        "mutations": [
            (30,'R'), (65,'T'), (72,'A'), (80,'R'), (143,'G'),
            (147,'P'), (157,'G'), (163,'A'), (167,'T'), (171,'V')
        ],
        "rationale": "S30R for folding + brightness mutations"
    },
    "Seq3_sfGFP_BrightPlus": {
        "backbone": sfGFP,
        "desc": "sfGFP + brightness enhancers",
        "mutations": [
            (72,'A'), (143,'G'), (147,'P'), (157,'G'), (167,'T'), (171,'V')
        ],
        "rationale": "Superfolder backbone + extra brightness"
    },
    "Seq4_sfGFP_Stable": {
        "backbone": sfGFP,
        "desc": "Thermostability focus",
        "mutations": [
            (18,'E'), (147,'P'), (157,'G'), (167,'T'), (171,'V')
        ],
        "rationale": "Minimal changes for max stability"
    },
    "Seq5_avGFP_Conservative": {
        "backbone": avGFP,
        "desc": "Winner-proven only",
        "mutations": [
            (30,'R'), (64,'F'), (65,'T'), (72,'A'), (80,'R'),
            (147,'P'), (163,'A'), (167,'T'), (171,'V'), (202,'D')
        ],
        "rationale": "Only mutations found in ≥2 previous winning sequences"
    },
    "Seq6_avGFP_Evolved": {
        "backbone": avGFP,
        "desc": "Diverse combination",
        "mutations": [
            (65,'T'), (72,'A'), (80,'R'), (143,'G'),
            (157,'G'), (167,'T'), (171,'V'), (203,'I')
        ],
        "rationale": "Includes T203I for spectral diversity"
    }
}

if __name__ == "__main__":
    print("="*60)
    print("FoldSynth Sequence Design Pipeline")
    print("="*60)
    
    results = []
    for name, info in DESIGN_STRATEGIES.items():
        seq = mutate(info["backbone"], info["mutations"])
        if not validate(seq):
            print(f"❌ {name}: Validation failed!")
            continue
        
        muts = []
        ref = info["backbone"]
        for i in range(len(ref)):
            if ref[i] != seq[i]:
                muts.append(f"{ref[i]}{i+1}{seq[i]}")
        
        print(f"\n✅ {name}")
        print(f"  Strategy: {info['desc']}")
        print(f"  Backbone: {'avGFP' if info['backbone']==avGFP else 'sfGFP'}")
        print(f"  Length: {len(seq)}")
        print(f"  Mutations ({len(muts)}): {', '.join(muts)}")
        results.append((name, seq))

    # Generate submission CSV
    with open('../submission.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Team_Name', 'Seq_ID', 'Sequence'])
        for i, (name, seq) in enumerate(results):
            writer.writerow(['FoldSynth', f'Seq{i+1}', seq])
    
    print(f"\n✅ Generated submission.csv with {len(results)} sequences")
