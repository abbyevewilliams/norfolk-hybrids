#!/bin/bash
#SBATCH --time=06:00:00               
#SBATCH --job-name=stats_zalb
#SBATCH --partition=short

#########################################################################################################
# SLIDING WINDOW STATISTICS USING THE GENOMICS_GENERAL REPO
# A script to convert VCF to geno format, then run ABBA BABA statistics in sliding windows.
# We can also run general pop gen stats.
#
# Abby Williams, November 2025
# abigail.williams@biology.ox.ac.uk
#########################################################################################################

# Set params
WIN_SIZE=50000
OVERLAP=10000
MIN_SITES=50
INPUT_VCF=/data/biol-silvereye/ball6625/norfolk-analyses/hard_call/results/hard_call_min_prob_0.80_min_percent_90.reheader.variantIDs.sorted.vcf.gz
BASE=$(basename $INPUT_VCF ".reheader.variantIDs.sorted.vcf.gz")
GENO=${BASE}.geno
OUT_DIR="./results"
OUTPUT_BASENAME="${OUT_DIR}/unpruned_hard_call_min_prob_0.80_min_percent_90_${WIN_SIZE}_${OVERLAP}_${MIN_SITES}"
mkdir -p $OUT_DIR

# Create input .geno file if it doesn't already exist
if [ ! -f $GENO]; then
	python3 ./genomics_general/VCF_processing/parseVCF.py -i $INPUT_VCF -o $GENO
fi

# ABBA-BABA (zten)
# Note that we used the -f phased option here for compatibility even though our data isn't phased
python ./genomics_general/ABBABABAwindows.py -w $WIN_SIZE -s $OVERLAP -m $MIN_SITES \
	-g $GENO -o "${OUTPUT_BASENAME}.zten.ABBABABA.csv" \
	--popsFile pops.txt -f phased \
	-P1 nz_zlat -P2 ni_zlat -P3 zten_hist -O Outgroup 

#ABBA-BABA (zalb)
python ./genomics_general/ABBABABAwindows.py -w $WIN_SIZE -s $OVERLAP -m $MIN_SITES \
	-g $GENO -o "${OUTPUT_BASENAME}.zalb.ABBABABA.csv" \
	--popsFile pops.txt -f phased \
	-P1 nz_zlat -P2 ni_zlat -P3 zalb -O Outgroup

# General pop gen statistics
#python ./genomics_general/popgenWindows.py -w $WIN_SIZE -s $OVERLAP -m $MIN_SITES \
#	-g $GENO -o "${OUTPUT_BASENAME}.stats.csv" \
# 	--popsFile pops.txt \
#	-p ni_zlat -p nz_zlat -p zten -p zalb \
#	-p zlat_hist -p zten_hist -p zlat_zten \
#	-p zlat_zalb -p Outgroup