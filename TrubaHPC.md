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
 - You need to have conda for package manager and snakemake for workflow batch submission in your HPC cluster. The following set of commands are given for installation conda in your environment without super user right. It might be useful if you need a specific version of snakemake or any required tool as well. Please have a look at the [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more information about Python3 version of Miniconda, mamba and snakemake installation.

**to install miniconda in your home folder** \
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh . \
$ bash Miniconda3-latest-Linux-x86_64.sh  \
$ source ~/miniconda3/etc/profile.d/conda.sh \
$ conda --help 

**mamba is officially recommended to install snakemake. We need to have mamba first.**\
$ conda install -c conda-forge mamba 

**to install snakemake via mamba**\
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake=5.32.0 

**Conda and snakemake are ready for job submission. Run followings to activate snakemake for each login**\
$ source ~/miniconda3/etc/profile.d/conda.sh \
$ conda activate snakemake \
$ snakemake --version 
 
**Caching between-workflow is an important functionality which avoids redundant computations if rules have already performed**\
- Snakemake 5.32 seems to not have a support caching between-workflow for cluster submission. Therefore, adding --cache parameter at 666th line on the following python file is required. Additionally, the path, hashed of input-sw-parameters are stored, should be set.\
$ vi ~/miniconda3/envs/snakemake/lib/python3.9/site-packages/snakemake/executors/\_\_init\_\_.py\
                    "--nocolor **--cache** --notemp --no-hooks --nolock "
- $ export SNAKEMAKE_OUTPUT_CACHE=/truba/home/emrah/shared/snakemake-cached
- Finally, getting previously computed output using caching on cluster, you need to make following change (change if not -> if) at 798th line on ~/miniconda3/envs/snakemake/lib/python3.9/site-packages/snakemake/io.py file. Otherwise, it raises MissingOutputException because of symbolic link.\
		    "if (is_flagged(f, "pipe") and ignore_pipe)"   

 - If you submit batchs for levrek1, Truba HPC, you have to do following steps, for blastb and paml. For other HPC, you need to put required files manually.\

## Blastdb
- You need to put **all_eu.fasta** file under **resources/blastdb** folder for alignment. \
$ cd resources \
$ rm blastdb \
$ ln -s /truba/home/emrah/shared/blastdb blastdb \
*These are the default path for blastdb and if you want, you can change both folder and its name in the content of config/config.yml and put the dbfile whereever you want.*

## Paml
- Although conda package manager installs all required software, you have to compile codeml manually under the path resources/paml4.9j. There is a guide how to install codeml that can be reached by [following this url](http://abacus.gene.ucl.ac.uk/software/paml.html).

`resources/paml4.9j should be accessible.` 

**NOTES:** There is already a compiled version of PAML (paml4.9j) for levrek1 (TRUBA HPC). You can easily use the installation (`/truba/home/emrah/shared/paml4.9j`) by running following commands. \
$ rm paml4.9j \
$ ln -s /truba/home/emrah/shared/paml4.9j paml4.9j

# Configuration
## Parameters
The content of config.yml under config folder indicates the name of proteins analayzed and parameters for all consequtive tasks performed.
- **workdir** indicates the working directory. After cloning the repository, you need to set the path (PWD) properly. Default is /cta/users/eakkoyun/WORKFOLDER/temp/phylogeny-workflow
- query_fasta lists all proteins that will be analyzed. All msa files should be stored resources/query_fasta. A few example of msa files are available in the repository.
- All other parameters for variety number of rules inside Snakefile. You can easily changes the parameter.
## Cluster
- There is a single file (config/slurm_truba/config.yaml) for batch submission. You need to make proper changes for your HPC environment. An example for Truba HPC (config.yaml) is given, please use proper partition name for your requirements.
```
jobs: 64 
cluster: "sbatch -p mid2 -J {rule}.job --qos long_investor --nodes=1 --ntasks=1 -t {resources.time_min}  -c {resources.cpus} -o logs/cluster/{rule}_%A.out -e logs/cluster/{rule}_%A.err"
default-resources: [cpus=4, time_min=9600]
```
# Run
- Running workflow on a HPC cluster is simply now. First, run snakemake with dry-run parameter to check that everything is fine. Then, delete the parameter and run the snakemake as following. It will submit job per rule defined in the Snakefile.

**$ cd phylogeny-snakemake** \
**$ pwd** # set workdir inside config/config.yml file with this path \
**$ cd workflow** \
**$ snakemake --use-conda --cache --profile ../config/slurm_truba --dry-run**


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
**$ snakemake --use-conda --cache --profile ../config/slurm_truba --keep-going --wms-monitor http://ephesus.sabanciuniv.edu:5000** \
Please pay attention to the following points for running snakemake workflow on HPC.
- use keep-going parameters to proceed running independent jobs in case of any failure on a task. This option allows submitting jobs for other proteins, while the consequtive jobs are not submitted for the protein that we observed a failure.
- use screen or execute the snakemake command in background. Otherwise, the next jobs are not submitted when we close the terminal or lost connection to the user interface. This is especially useful for long set of runs in workflow. If you do not know how to use screen in linux, you can execute the command as following and follow the output under .snakemake/log/ folder.\
- use panoptes server for monitoring submitted jobs. The snakemake std outputs during execution is transformed into gui which can be visualized at the given url.

**$ nohup snakemake --use-conda --cache --profile ../config/slurm_truba --keep-going --wms-monitor http://ephesus.sabanciuniv.edu:5000 > /dev/null 2>&1 &**

# Run multiple protein within single node (optional)
- Truba has a new cluster (hamsi), where users submit their jobs requesting 28 cores so we have to run creating ml_tree for multiple proteins within a job.
- Submitting jobs requesting only 4 cores are allowed to run mid2 partition. The other tasks in the workflow does not take huge amount of CPU computation time. They can be run on mid2 partitions until ml_tree as following.\
  **$ snakemake --use-conda --cache --profile ../config/slurm_truba --keep-going --until trim_msa**
- Jobs creating ml tree should be submitted to hamsi partition. Each ml_tree task request 7 cpus, it means that 4 protein can be joined within the single job as following. Do not forget to change partition name to hamsi.\
**$ snakemake --use-conda --cache --profile ../config/slurm_truba --keep-going --groups somerule=ml_tree --group-components ml_tree=4 --until ml_tree**








