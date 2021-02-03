# Phylogeny workflow

This repository prepared for batch submission to an HPC cluster (slurm) and running following workflow. You can access the DAG file by [clicking here](https://raw.githubusercontent.com/emrahkyn/phylogeny-workflow/main/images/rulegraph.svg?token=ASHDZ5HNMJ73KOJ23CFSQ3K773LUE).
1. blastp
2. get_blasthits
3. header_update
4. msa
5. trim_msa
6. remove_gaps
7. ml_tree
8. unroot_tree
9. run_codeml
10. compute_score

The workflow is mainly prepared for tosun [@Sabanci HPC](https://github.com/CompGenomeLab/phylogeny-snakemake/blob/main/SabanciHPC.md) and levrek1 [@Truba HPC](https://github.com/CompGenomeLab/phylogeny-snakemake/blob/main/TrubaHPC.md). However, you can also use the repository for any other HPC, cloud computing platform by updating the config file (profile).
