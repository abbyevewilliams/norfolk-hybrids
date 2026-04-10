#!/usr/bin/env python3
"""
strip_invariants_iqtree_compatible.py
Keep only columns that IQ-TREE would consider variable.

Rule implemented:
 - For each column, map characters to IUPAC base sets (A,C,G,T, plus ambiguity codes).
 - Remove gap characters '-' and '.' from consideration.
 - If the intersection of all character-sets (over non-gap characters) is non-empty,
   then the column CAN be explained by a single nucleotide -> treat as INVARIANT (remove).
 - Otherwise keep the column (variable).
Outputs diagnostics to stdout (counts + first few invariant column descriptions).
"""
import sys
from Bio import SeqIO

IUPAC_SETS = {
    'A': set('A'),
    'C': set('C'),
    'G': set('G'),
    'T': set('T'),
    'U': set('T'),
    'R': set('AG'),
    'Y': set('CT'),
    'S': set('GC'),
    'W': set('AT'),
    'K': set('GT'),
    'M': set('AC'),
    'B': set('CGT'),
    'D': set('AGT'),
    'H': set('ACT'),
    'V': set('ACG'),
    'N': set('ACGT'),
    # lower-case support
    'a': set('A'),
    'c': set('C'),
    'g': set('G'),
    't': set('T'),
    'u': set('T'),
    'r': set('AG'),
    'y': set('CT'),
    's': set('GC'),
    'w': set('AT'),
    'k': set('GT'),
    'm': set('AC'),
    'b': set('CGT'),
    'd': set('AGT'),
    'h': set('ACT'),
    'v': set('ACG'),
    'n': set('ACGT'),
    '.': set(),  # treat as gap (ignored)
    '-': set(),  # gap
}

def char_to_set(ch):
    if ch in IUPAC_SETS:
        return IUPAC_SETS[ch]
    # fallback: if unexpected char, treat as N (conservative)
    return set('ACGT')

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python strip_invariants_iqtree_compatible.py <input.fasta> <output.fasta>")

    infile, outfile = sys.argv[1], sys.argv[2]
    recs = [r for r in SeqIO.parse(infile, "fasta")]
    if not recs:
        sys.exit("ERROR: no sequences read from input")

    nseq = len(recs)
    L = len(recs[0].seq)
    for r in recs[1:]:
        if len(r.seq) != L:
            sys.exit("ERROR: sequences are not equal length")

    mask = [False] * L
    invariant_cols = []
    kept = 0

    for i in range(L):
        col_chars = [str(r.seq)[i] for r in recs]
        # build list of sets for non-gap characters
        sets = []
        for ch in col_chars:
            if ch in ('-', '.'):
                continue
            s = char_to_set(ch)
            sets.append(s)
        # if no non-gap characters -> treat as invariant (all gaps)
        if not sets:
            invariant_cols.append((i, "<all-gaps>"))
            continue
        # compute intersection
        inter = set.intersection(*sets)
        if inter:
            # there is at least one base compatible with every symbol -> invariant
            desc = "/".join(sorted({c for ch in set(col_chars) if ch not in ('-','.') for c in char_to_set(ch)}))
            invariant_cols.append((i, desc))
            continue
        else:
            mask[i] = True
            kept += 1

    # write output
    with open(outfile, "w") as out:
        for r in recs:
            seq_out = ''.join(str(r.seq)[i] for i,keep in enumerate(mask) if keep)
            out.write(f">{r.id}\n{seq_out}\n")

    total = L
    inv_count = total - kept
    print(f"Input columns: {total}")
    print(f"Kept (variable) columns: {kept}")
    print(f"Invariant columns (per IQ-TREE compatibility rule): {inv_count}")
    print(f"Number of sequences: {nseq}")

if __name__ == '__main__':
    main()
