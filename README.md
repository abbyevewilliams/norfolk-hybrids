# Code for reproducing analyses from *Recent hybridisation and ghost introgression among a trio of island passerines.*

Abby Williams, April 2026

Email: abigail.williams (at) biology.ox.ac.uk

In this repo:
* popglen_config.yaml: config file fpr running [PopGLen](https://github.com/zjnolen/PopGLen)
* hard_call/: code for producing a 'hard call' dataset from ANGSD genotype likelihoods (beagle file).
* plink_commands.txt: list of plink commands used on 'hard call' dataset.
* phylo/: code for reproducing phylogenomic analyses.
* hyde.sh: shell script for running [HyDe](https://github.com/pblischak/HyDe/tree/main) (hybridisation analyses).
* dsuite.sh: shell script for running [Dsuite](https://github.com/millanek/Dsuite) (introgression analyses).
* stats.sh: shell script for running [genomics_general](https://github.com/simonhmartin/genomics_general)
* extract_genes.sh: shell script for extracting genes from top 1% of introgressed windows and extracting gene candidates in those windows.
* triangular.Rmd: R markdown file for running TriangulaR package (https://github.com/omys-omics/triangulaR)

