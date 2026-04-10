# ---------------------------------------------------
#
# A script to convert beagle to VCF whilst filtering
# for sites where more than X% of individuals have a 
# GP of >Y%.
#
# Based on a script by Andrea Estandia 
# (github.com/andreaestandia)
#
# Usage: python beagle2VCF.py input.beagle output.vcf \
# --min-prob 0.9 --min-percent 90 \
# --fai path/to/fai
#
# ---------------------------------------------------

import sys
from datetime import datetime
import argparse
import gzip

def translate_alleles(allele1, allele2):
    mapping = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
    try:
        return mapping[allele1], mapping[allele2]
    except KeyError:
        raise ValueError(f"Unrecognized allele: {allele1} or {allele2}")

def parse_args():
    parser = argparse.ArgumentParser(description="Convert BEAGLE to filtered VCF.")
    parser.add_argument("input", help="Input BEAGLE file")
    parser.add_argument("output", help="Output VCF file")
    parser.add_argument("--min-prob", type=float, default=0.9,
                        help="Minimum genotype probability (default: 0.9)")
    parser.add_argument("--min-percent", type=float, default=100,
                        help="Minimum percent of samples that must meet min-prob (default: 100)")
    parser.add_argument("--fai", help="Reference FASTA index (.fai) file with contig lengths", 
                        required=True)
    return parser.parse_args()

def open_file(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding='utf-8')
    else:
        return open(filename, mode, encoding='utf-8' if 't' in mode else None)

def read_fai(fai_path):
    contigs = {}
    with open(fai_path, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            contigs[fields[0]] = int(fields[1])
    return contigs

args = parse_args()

sample_names = []
total_sites = 0
written_sites = 0

now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
sys.stderr.write(f"START TIME: [{now}]\n")
sys.stderr.write(f"\nMinimum genotype probability: {args.min_prob}\n")
sys.stderr.write(f"Minimum percentage of individuals required: {args.min_percent}\n")
sys.stderr.write(f"Input beagle file: {args.input}\n")
sys.stderr.write(f"Input .fai: {args.fai}\n")

# Read .fai
contig_lengths = read_fai(args.fai)

with open_file(args.input, 'rt') as infile, open_file(args.output, 'wt') as outfile:
    # Read header
    for line in infile:
        if line.startswith('marker'):
            sample_names = line.strip().split()[3::3]
            break

    if not sample_names:
        sys.stderr.write("ERROR: Could not find sample names in BEAGLE file.\n")
        sys.exit(1)

    n_samples = len(sample_names)

    # Write VCF header
    outfile.write('##fileformat=VCFv4.2\n')
    outfile.write(f'##min_prob={args.min_prob}\n')
    outfile.write(f'##min_percent={args.min_percent}\n')
    outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

    # Print contig lines sorted by name or order in fai
    for contig, length in contig_lengths.items():
        outfile.write(f'##contig=<ID={contig},length={length}>\n')

    outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    for name in sample_names:
        outfile.write(f"\t{name}")

    # Process each line
    line_num = 0
    for line in infile:
        line_num += 1
        total_sites += 1
        cols = line.strip().split()

        if len(cols) != (3 + 3 * n_samples):
            #sys.stderr.write(f"WARNING: Line {line_num} skipped due to column mismatch.\n")
            continue

        try:
            ref, alt = translate_alleles(cols[1], cols[2])
            chrom, pos = cols[0].split('_')
        except Exception as e:
            sys.stderr.write(f"WARNING: Line {line_num} skipped: {e}\n")
            continue

        genotypes = []
        n_passing = 0

        try:
            for i in range(3, len(cols), 3):
                probs = list(map(float, cols[i:i+3]))
                max_prob = max(probs)
                if max_prob >= args.min_prob:
                    n_passing += 1

                if probs[0] > probs[1] and probs[0] > probs[2]:
                    gt = '0/0'
                elif probs[1] > probs[0] and probs[1] > probs[2]:
                    gt = '0/1'
                else:
                    gt = '1/1'

                genotypes.append(gt)
        except Exception as e:
            sys.stderr.write(f"WARNING: Line {line_num} skipped due to GP parse error: {e}\n")
            continue

        percent_passing = (n_passing / n_samples) * 100
        if percent_passing >= args.min_percent:
            outfile.write(f"\n{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT")
            for gt in genotypes:
                outfile.write(f"\t{gt}")
            written_sites += 1
        #else:
            #sys.stderr.write(f"INFO: Line {line_num} skipped ({percent_passing:.1f}% ≥ {args.min_prob}).\n")

# Final summary
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
sys.stderr.write(f"\nEND TIME: [{now}]\n")
sys.stderr.write(f"\nSUMMARY:\n")
sys.stderr.write(f"  Total sites processed: {total_sites}\n")
sys.stderr.write(f"  Sites written to VCF:   {written_sites}\n")
sys.stderr.write(f"  Sites filtered out:     {total_sites - written_sites}\n")
