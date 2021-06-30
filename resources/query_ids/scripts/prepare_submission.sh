#source prepare_submission.sh sabanci

cluster_name=$1
if [ $cluster_name == 'sabanci' ]; then
	python3 prepare_config.py sabanci > ../../../config/config.yml 2>/dev/null
	export SNAKEMAKE_OUTPUT_CACHE=/cta/groups/adebali/static/snakemake-cached
	conda activate snakemake
	echo "Config is ready for submission. To run snakemake on $cluster_name HPC: " 
	echo "cd ../../../workflow && nohup snakemake --use-conda --cache --profile ../config/slurm_sabanci --keep-going --wms-monitor http://ephesus.sabanciuniv.edu:5000 > /dev/null 2>&1 &"
elif [ $cluster_name == 'truba' ]; then
	python3 prepare_config.py truba > ../../../config/config.yml 2>/dev/null
	rm -rf ../../blastdb && ln -s /truba/home/emrah/shared/blastdb ../../blastdb
	rm -rf ../../paml4.9j && ln -s /truba/home/emrah/shared/paml4.9j ../../paml4.9j
	export SNAKEMAKE_OUTPUT_CACHE=/truba/home/emrah/shared/snakemake-cached
	source ~/miniconda3/etc/profile.d/conda.sh
	conda activate snakemake
	echo "Config is ready for submission. To run snakemake on $cluster_name HPC: "
	echo "cd ../../../workflow && nohup snakemake --use-conda --cache --profile ../config/slurm_truba --keep-going > /dev/null 2>&1 &" 
fi
