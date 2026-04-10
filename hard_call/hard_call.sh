#!/bin/bash
#SBATCH --time=01:00:00            
#SBATCH --job-name=hard_call     
#SBATCH --partition=short
#SBATCH --output=sbatch_%j.log
#SBATCH --error=sbatch_%j.error

#########################################################################################################
# CONVERT BEAGLE FILES TO "HARD CALL" VCF
# The input for this script is the .beagle.gz genotype likelihood files outputted by the PopGLen pipeline.
# First we combine all the .beagle.gz "chunks" into a single .beagle.gz file.
# Then we use a custom python script, beagle2VCF.py, to filter out low-confidence genotype likelihoods.
# This script saves the output to VCF, taking the high-confidence genotype likelihoods as hard calls.
# Finally, we reheader the VCF using a .txt file of sample names, and annotate variant IDs for compatibility with downstream tools.

# Abby Williams, August 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

echo `date`
set -euo pipefail

# ---- Load software (BCFtools v1.14) ----
ml BCFtools/1.14-GCC-11.2.0

# ---- Set input paths ----
# CHUNK_DIR: Directory containing beagle chunks
# INPUT_BEAGLE: Name of combined beagle we will input to beagle2VCF.py
# REF_INDEX: .fai file for hard calling
# SAMPLE_NAMES: tab separated file containing old sample names (col 1) and new sample names (col 2)
CHUNK_DIR=/data/biol-silvereye/ball6625/popglen/results/datasets/main/beagles/chunks
INPUT_BEAGLE=./combined.beagle.gz
REF_INDEX=/data/biol-silvereye/ball6625/ref-genome/zosterops_lateralis/GCA_965231275.1_bZosLat1.hap1.1_genomic.fna.fai
SAMPLE_NAMES=./reheader.txt
echo "CHUNK_DIR: ${CHUNK_DIR} INPUT_BEAGLE: ${INPUT_BEAGLE}"

# ---- Set params ----
# MIN_PROB: the minimum posterior probability in a particular genotype likelihood
# MIN_PERCENT: the minimum percentage of individuals needing to have MIN_PROB or above
MIN_PROB=0.80
MIN_PERCENT=90
HARD_CALL_VCF=hard_call_min_prob_${MIN_PROB}_min_percent_${MIN_PERCENT}.vcf
echo "MIN_PROB: ${MIN_PROB} MIN_PERCENT ${MIN_PERCENT} HARD_CALL_VCF: ${HARD_CALL_VCF}"

# ---- Concatenate all beagle files ----
if [ ! -f $INPUT_BEAGLE ]; then
    echo "Combining beagle files..."
    (zcat $CHUNK_DIR/*.beagle.gz \
     | grep -m 1 "^marker"; zcat $CHUNK_DIR/*.beagle.gz \
     | grep -v "^marker") | bgzip > $INPUT_BEAGLE
fi

# ---- Run hard calling ----
if [ ! -f $HARD_CALL_VCF ]; then
  echo "Running hard calls..."
  python ./beagle2VCF.py ${INPUT_BEAGLE} ${HARD_CALL_VCF} \
      --min-prob ${MIN_PROB} --min-percent ${MIN_PERCENT} \
      --fai ${REF_INDEX}
  bgzip ${HARD_CALL_VCF}
fi

# ---- Reheader VCF ----
echo "Reheadering..."
VCF_PREFIX="$(basename "$HARD_CALL_VCF" .vcf)"
# use -Ou (uncompressed BCF) when piping for speed / lower mem
bcftools view -Ou "${HARD_CALL_VCF}.gz" \
  | bcftools reheader -s "$SAMPLE_NAMES" \
  | bcftools view -Oz -o "${VCF_PREFIX}.reheader.vcf.gz"
tabix -p vcf "${VCF_PREFIX}.reheader.vcf.gz"

# ---- Annotate variant IDs (only sets ID where missing) ----
echo "Annotating IDs..."
zcat "${VCF_PREFIX}.reheader.vcf.gz" \
  | awk 'BEGIN{OFS="\t"} /^#/ {print; next}
         {
           # take first ALT if there are multiple
           split($5, a, ",");
           alt = a[1];
           if ($3 == ".") $3 = $1"_"$2;
           print
         }' \
  | bgzip -c > "${VCF_PREFIX}.reheader.variantIDs.vcf.gz"

# ---- Sort the annotated VCF and write output ----
echo "Sorting annotated VCF..."
bcftools sort -Oz -o "${VCF_PREFIX}.reheader.variantIDs.sorted.vcf.gz" \
  "${VCF_PREFIX}.reheader.variantIDs.vcf.gz"
tabix -p vcf "${VCF_PREFIX}.reheader.variantIDs.sorted.vcf.gz"

# ---- Remove intermediates ----
echo "Removing intermediates..."
rm ${INPUT_BEAGLE} ${HARD_CALL_VCF}
echo "Done"