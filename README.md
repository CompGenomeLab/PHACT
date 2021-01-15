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

# Pre-requirements
## Snakemake and Conda
 - You need to have conda for package manager and snakemake for workflow batch submission in your HPC cluster. 
## Blastdb
- You need to put **all_eu.fasta** file under **resources/blastdb** folder for alignment. 
*These are the default path for blastdb and if you want, you can change both folder and its name in the content of config/config.yml and put the dbfile whereever you want.*

## Paml
- Although conda package manager installs all required software, you have to compile codeml manually under the path resources/paml4.9j. There is a guide how to install codeml that can be reached by [following this url](http://abacus.gene.ucl.ac.uk/software/paml.html).

`resources/paml4.9j/bin/codeml.exe should be accessible.` 

**NOTES:** There is already a compiled version of PAML (paml4.9j) for tosun (Sabanci HPC). You can easily copy the folder (`/cta/users/eakkoyun/WORKFOLDER/paml4.9j`) into your environment.

# Configuration
## Parameters
The content of config.yml under config folder indicates the name of proteins analayzed and parameters for all consequtive tasks performed.
- output_folder indicates the working directory. After cloning the repository, you need to set the path (PWD) properly. Default is /cta/users/eakkoyun/WORKFOLDER/temp/phylogeny-workflow
- msa_name lists all proteins that will be analyzed. All msa files should be stored resources/msa_files. A few example of msa files are available in the repository.
- All other parameters for variety number of rules inside Snakefile. You can easily changes the parameter.
## Cluster
- There is a single file (config/slurm/config.yaml) for batch submission. You need to make proper changes for your HPC environment. You do not have to make changes for Sabanci HPC cluster.







