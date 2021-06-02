#!/bin/bash

trimmed_msa=$1
raxml_model=$2
raxml_out_name=$3
raxml_seed=$4
output_tree=$5
tree_num=$6

# recommended # threads
thread_num=`raxml-ng --parse --msa $trimmed_msa --model $raxml_model --nofile|grep "Recommended number of threads"|awk -F: {'print $2'}| xargs`
echo "Recommended number of thread: $thread_num"
echo "Number of CPUs on the allocated nodes: $SLURM_CPUS_ON_NODE" 

if (("$thread_num" <= "$SLURM_CPUS_ON_NODE")); then
	folder=`dirname $output_tree`
	raxml-ng  --msa $trimmed_msa --data-type AA --model $raxml_model --prefix $folder/$raxml_out_name --seed $raxml_seed --threads $thread_num --tree rand{$tree_num} pars{$tree_num}
else
	echo "ERROR: Required CPU number: $thread_num. Assign more CPUs for the next submission"
fi

