#./archive_completed_query.sh truba

finished_num_task=5
cluster_name=$1
if [ $cluster_name == 'sabanci' ]; then
	archived_path=/cta/groups/adebali/static/archived
elif [ $cluster_name == 'truba' ]; then
	archived_path=/truba/home/emrah/shared/archived
fi

count_completed=0
count_incompleted=0
for protein in `ls ../../../results/`; do 
	num_tasks=`ls ../../../results/$protein|wc -l`; 
        if [ "$num_tasks" -eq "$finished_num_task" ]; then
		if [ ! -d "$archived_path/$protein" ]; then 
			cp -R -L ../../../results/$protein $archived_path/$protein
			echo "$protein archived"
		fi
                let "count_completed=count_completed+1"
	else
		let "count_incompleted=count_incompleted+1"
	fi
done
echo "Completed protein num: $count_completed"
echo "Not computed yet protein num: $count_incompleted"
