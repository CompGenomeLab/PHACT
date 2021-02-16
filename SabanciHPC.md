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
 - You need to clone the repository in your working directory.
 `git clone https://github.com/CompGenomeLab/phylogeny-snakemake.git`
 - You need to have conda for package manager and snakemake for workflow batch submission in your HPC cluster. Loading snakemake module is enough for Sabanci HPC. 
 `module load snakemake-5.23.0`
 - If you need to install snakemake into your environment, for example, latest version of snakemake, please have a look at the [installation link](https://github.com/CompGenomeLab/phylogeny-snakemake/blob/main/TrubaHPC.md#snakemake-and-conda). Since the installed version of snakemake does not support caching between workflow, you must follow the installation instructions to deploy snakemake in your home folder, and make a simple change on the official python code. Please set the cached path as following for Sabanci HPC \
 `export SNAKEMAKE_OUTPUT_CACHE=/cta/groups/adebali/static/snakemake-cached`
 - If you submit batchs for tosun, Sabanci HPC, you do not have to do following steps, for blastb and paml.
## Blastdb
- You need to put **all_eu.fasta** file under **resources/blastdb** folder for alignment. 
*These are the default path for blastdb and if you want, you can change both folder and its name in the content of config/config.yml and put the dbfile whereever you want.*

## Paml
- Although conda package manager installs all required software, you have to compile codeml manually under the path resources/paml4.9j. There is a guide how to install codeml that can be reached by [following this url](http://abacus.gene.ucl.ac.uk/software/paml.html).

`resources/paml4.9j/bin/codeml.exe should be accessible.` 

**NOTES:** There is already a compiled version of PAML (paml4.9j) for tosun (Sabanci HPC) and symbolic linked for your usage.

# Configuration
## Parameters
The content of config.yml under config folder indicates the name of proteins analayzed and parameters for all consequtive tasks performed.
- **workdir** indicates the working directory. After cloning the repository, you need to set the path (PWD) properly. Default is /cta/users/eakkoyun/WORKFOLDER/temp/phylogeny-workflow
- query_fasta lists all proteins that will be analyzed. All msa files should be stored resources/query_fasta. A few example of msa files are available in the repository.
- All other parameters for variety number of rules inside Snakefile. You can easily changes the parameter.
## Cluster
- There is a single file (config/slurm/config.yaml) for batch submission. You need to make proper changes for your HPC environment. You do not have to make changes for Sabanci HPC cluster.
# Run
- Running workflow on a HPC cluster is simply now. First, run snakemake with dry-run parameter to check that everything is fine. Then, delete the parameter and run the snakemake as following. It will submit job per rule defined in the Snakefile.

**$ cd phylogeny-snakemake** \
**$ pwd** # set workdir inside config/config.yml file with this path \
**$ cd workflow** \
**$ snakemake --use-conda --cache --profile ../config/slurm_sabanci --dry-run**


```
Job counts:
	count	jobs
	1	all
	3	blastp
	3	compute_score
	3	get_blasthits
	3	header_update
	3	ml_tree
	3	msa
	3	remove_gaps
	3	run_codeml
	3	trim_msa
	3	unroot_tree
	31
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```
**$ snakemake --use-conda --cache --profile ../config/slurm_sabanci --keep-going --wms-monitor http://ephesus.sabanciuniv.edu:5000** \
Please pay attention to the following points for running snakemake workflow on HPC.
- use keep-going parameters to proceed running independent jobs in case of any failure on a task. This option allows submitting jobs for other proteins, while the consequtive jobs are not submitted for the protein that we observed a failure.
- use screen or execute the snakemake command in background. Otherwise, the next jobs are not submitted when we close the terminal or lost connection to the user interface. This is especially useful for long set of runs in workflow. If you do not know how to use screen in linux, you can execute the command as following and follow the output under .snakemake/log/ folder.\
- use panoptes server for monitoring submitted jobs. The snakemake std outputs during execution is transformed into gui which can be visualized at the given url.

**$ nohup snakemake --use-conda --cache --profile ../config/slurm_sabanci --keep-going  --wms-monitor http://ephesus.sabanciuniv.edu:5000 > /dev/null 2>&1 &**







