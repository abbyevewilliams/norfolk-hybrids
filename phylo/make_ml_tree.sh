#!/bin/bash
#SBATCH --time=2-00:00:00               
#SBATCH --job-name=run_iqtree
#SBATCH --partition=long
#SBATCH --output=%A_%a.log
#SBATCH --error=%A_%a.error
#SBATCH --mem-per-cpu=32G

#########################################################################################################
# MAXIMUM LIKELIHOOD TREE FROM SNP DATA
# A script to run IQTREE for all concatenated SNPs (thinned) to make the maximum likelihood tree.
#
# Abby Williams, November 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

# Set paths
IQ_OUTDIR=./results
mkdir -p $IQ_OUTDIR
VCF=/data/biol-silvereye/ball6625/norfolk-analyses/plink/04_1_thinned.vcf

# Params
OUTPUT_PREFIX=thinned
# MIN_SAMPLES is the minimum number of samples required for a site to be written in the .fasta file
MIN_SAMPLES=4

# Load software (IQTREE 2.4)
ml Anaconda3
source activate ./env

# Convert to fasta format
FASTA_OUTPUT=${OUTPUT_PREFIX}.min${MIN_SAMPLES}.fasta
if [ ! -f $FASTA_OUTPUT ]; then
	echo "Converting to fasta..."
	python vcf2phylip.py -i $VCF --output-prefix $OUTPUT_PREFIX -f --phylip-disable -m $MIN_SAMPLES
fi

# Remove invariant sites
echo "Removing invariant sites..."
python strip_invariant_sites.py $FASTA_OUTPUT $FASTA_OUTPUT.invar_removed

# Run IQTREE
echo "Running IQTREE..."
iqtree -s $FASTA_OUTPUT.invar_removed -st DNA \
	-m MFP+ASC -bb 1000 -alrt 1000 -bnni -nt 1 \
	-o SRR12615514_SRR12615514 \
	-pre $OUTPUT_PREFIX