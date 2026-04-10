#!/bin/bash
#SBATCH --time=01:00:00               
#SBATCH --job-name=hyde
#SBATCH --partition=short
#SBATCH --output=%A_%a.log
#SBATCH --error=%A_%a.error

#########################################################################################################
# HYDE FOR HYBRID DETECTION
# A script to convert a VCF file to phylip, then run HyDe to detect hybrids.
#
# Abby Williams, November 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

# Load software
ml Anaconda3
ml BCFtools/1.14-GCC-11.2.0
source activate ./env

# Set paths
VCF=/data/biol-silvereye/ball6625/norfolk-analyses/plink/01_2_hwe_filtered.vcf

# Params
OUTPUT_PREFIX=all_samples
# MIN_SAMPLES is the minimum number of samples required for a site to be written in the phylip file
MIN_SAMPLES=4
N_IND=45
N_TAXA=7

# Convert to phylip format
PHYLIP_OUTPUT=${OUTPUT_PREFIX}.min${MIN_SAMPLES}.phy
if [ ! -f $PHYLIP_OUTPUT ]; then
	echo "Converting to phylip..."
	python vcf2phylip.py -i $VCF --output-prefix $OUTPUT_PREFIX -m $MIN_SAMPLES
fi

# Count number of sites in original VCF for HyDe
SITES=$(bcftools view -H $VCF | wc -l)

# Run HyDe (overall check)
run_hyde.py -i $PHYLIP_OUTPUT -m map.txt -o Outgroup -n $N_IND -t $N_TAXA -s $SITES

# Now check individuals
individual_hyde.py -i $PHYLIP_OUTPUT -m map.txt -tr triples.txt -n $N_IND -t $N_TAXA -s $SITES -o Outgroup