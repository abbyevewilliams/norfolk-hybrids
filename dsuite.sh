#!/bin/bash
#SBATCH --time=00:20:00               
#SBATCH --job-name=dsuite
#SBATCH --partition=short       
#SBATCH --output=dsuite_%A.log      
#SBATCH --error=dsuite_%A.error

#########################################################################################################
# RUN DSUITE
# The input for this script is the .vcf.gz "hard call" file.
# We run Dsuite on all possible combinations of populations to check for signatures of hybridisation.
# We also use the --ABBAclustering option to check whether significant windows cluster in the genome.

# Abby Williams, August 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

# Software setup
export PATH=/data/biol-silvereye/ball6625/software/Dsuite/Build:$PATH

# Set paths
INPUT_VCF=/data/biol-silvereye/ball6625/norfolk-analyses/hard_call/results/hard_call_min_prob_0.80_min_percent_90.reheader.variantIDs.sorted.vcf.gz
DATASET_NAME=$(basename $INPUT_VCF .reheader.variantIDs.sorted.vcf.gz)

# Run Dtrios
echo "Running Dsuite..."
Dsuite Dtrios \
	-n ${DATASET_NAME} \
	--ABBAclustering \
	${INPUT_VCF} \
	SETS.txt