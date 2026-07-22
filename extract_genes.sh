#!/bin/bash
#SBATCH --time=06:00:00               
#SBATCH --job-name=extract_genes_zalb
#SBATCH --partition=short
#SBATCH --output=extract_%A.log      
#SBATCH --error=extract_%A.error

#########################################################################################################
# EXTRACT GENE HITS FROM GENOMIC WINDOWS
# Maps windows from the Zosterops lateralis genome to annotated genes in the Taeniopygia guttata genome.
#
# Abby Williams, November 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

# Set paths
DATASET_NAME=zalb
INPUT_TSV=/data/biol-silvereye/ball6625/norfolk-analyses/gsea/input/zalb_top_1_percent_windows_50kb_10kb.tsv
REF=/data/biol-silvereye/ball6625/ref-genome/zosterops_lateralis/GCA_965231275.1_bZosLat1.hap1.1_genomic.fna
RESULTS_DIR=./results/${DATASET_NAME}
mkdir -p $RESULTS_DIR

# Load software
ml Anaconda3
source activate ./ncbi-env

# Convert input windows to fasta
echo "Converting input windows to fasta..."
tr -d '\r' < "$INPUT_TSV" \
  | awk 'BEGIN{OFS="\t"} NF>=3 {
      chrom=$1;
      start_raw=$2; end_raw=$3;
      # remove non-digits from coords but keep empty if missing
      gsub(/[^0-9]/,"",start_raw); gsub(/[^0-9]/,"",end_raw);
      start = (start_raw=="" ? 0 : start_raw-1); if(start<0) start=0;
      end   = (end_raw=="" ? start+1 : end_raw);
      if(end <= start) end = start + 1;
      print chrom, start, end, "win_" NR
    }' \
  | sort -k1,1 -k2,2n | uniq > "$RESULTS_DIR/windows.bed"
bedtools getfasta -fi "$REF" -bed "$RESULTS_DIR/windows.bed" -fo "$RESULTS_DIR/windows.fa" -name

# Mask repetitive regions in fasta
echo "Masking repetitive regions..."
dustmasker -in $RESULTS_DIR/windows.fa -out $RESULTS_DIR/windows.masked.fa -infmt fasta -outfmt fasta

# Query input fasta against zebra finch database
echo "Doing BLAST query..."
blastn -query "$RESULTS_DIR/windows.masked.fa" -db zebrafinch_db \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -evalue 1e-6 -max_target_seqs 50 \
  -out "$RESULTS_DIR/hits_raw.tsv"

# Filter and sort output based on identity, length of alignment, coverage of query
echo "Filtering BLAST output..."
awk 'BEGIN{FS=OFS="\t"}
  {
    pident=$3; alen=$4;
    if(pident>=60 && alen>=100) print $0
  }' "$RESULTS_DIR/hits_raw.tsv" > "$RESULTS_DIR/hits_filtered.tsv"

sort -k1,1 -k12,12nr "${RESULTS_DIR}/hits_filtered.tsv" > "${RESULTS_DIR}/hits_filtered_sorted.tsv"

# Only keep top 10 entries per query
echo "Selecting best entries..."
awk 'BEGIN{FS=OFS="\t"} {cnt[$1]++; if(cnt[$1]<=10) print $0}' "${RESULTS_DIR}/hits_filtered_sorted.tsv" > "${RESULTS_DIR}/hits_top10.tsv"

# Create a 5-column subject BED: chrom, start0, end, qid, bitscore
echo "Building BED file..."
awk 'BEGIN{FS=OFS="\t"}
{
  sid = $2;
  # strip common wrappers like ref|...| or gi|
  gsub(/^ref\|/,"",sid); sub(/\|$/,"",sid); gsub(/^gi\|/,"",sid);
  # get subject coords (sstart/send)
  s = ($9 < $10 ? $9 : $10);
  e = ($9 < $10 ? $10 : $9);
  if(s == "" || e == "" || s+0 == 0 && e+0 == 0) next;   # skip malformed lines
  if(s < 1) s = 1;
  start0 = s - 1;
  if(e <= start0) e = start0 + 1;
  qid = $1;
  bits = $12 + 0;    # bitscore (force numeric)
  print sid, start0, e, qid, bits;
}' "${RESULTS_DIR}/hits_top10.tsv" > "${RESULTS_DIR}/hits_top10.with_scores.bed"

# sort and uniq just in case
sort -k1,1 -k2,2n "${RESULTS_DIR}/hits_top10.with_scores.bed" | uniq > "${RESULTS_DIR}/hits_top10.with_scores.sorted.bed"

# Merge nearby subject hits (within 1kb), aggregating distinct query ids and mean bitscore
echo "Merging hits and aggregating stats..."
bedtools merge -i "${RESULTS_DIR}/hits_top10.with_scores.sorted.bed" -d 1000 -c 4,5 -o distinct,mean \
  > "${RESULTS_DIR}/hits_top10_merged.bed"

# Get gene names from zebra finch gff
echo "Getting gene names..."
bedtools intersect \
  -a "${RESULTS_DIR}/hits_top10_merged.bed" \
  -b zebra_finch_genes.gff \
  -wa -wb > "${RESULTS_DIR}/hits_with_genes.tsv"

# Extract gene names
awk -F'\t' '{
  attr = $NF;
  gene = "NA";
  if (match(attr,/Name=([^;]+)/,m)) gene=m[1];
  else if (match(attr,/gene_name=([^;]+)/,m)) gene=m[1];
  else if (match(attr,/gene=([^;]+)/,m)) gene=m[1];
  else if (match(attr,/ID=([^;]+)/,m)) gene=m[1];
  print $4,gene,$5;
}' OFS="\t" "${RESULTS_DIR}/hits_with_genes.tsv" > "${RESULTS_DIR}/hits_genes_named.tsv"

# Summarise
echo "Summarising..."
awk -F'\t' '$3 >= 500 {print}' "${RESULTS_DIR}/hits_genes_named.tsv" \
  | awk '!seen[$2]++' > "${RESULTS_DIR}/results.tsv"
awk '{print $2}' "${RESULTS_DIR}/results.tsv" > "${RESULTS_DIR}/results_short.txt"