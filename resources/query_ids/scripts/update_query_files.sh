#./update_query_files.sh truba  2>/dev/null 

cluster_name=$1
if [ $cluster_name == 'sabanci' ]; then
        archived_path=/cta/groups/adebali/static/archived
elif [ $cluster_name == 'truba' ]; then
        archived_path=/truba/home/emrah/shared/archived
fi


all_query_file=../all
completed_query_file=../completed
not_started_yet_query_file=../not_yet_computed
sabanci_query_file=../under_computation_tosun
truba_query_file=../under_computation_levrek1
finished_num_task=5

#empty the content of query files
cat /dev/null  > $completed_query_file
cat /dev/null  > $not_started_yet_query_file

count_completed=0
count_incompleted=0
count_under_computation=0
for protein in `cat $all_query_file`; do
        num_tasks=`ls $archived_path/$protein|wc -l`
        if [ "$num_tasks" -eq "$finished_num_task" ]; then
                echo "$protein" >> $completed_query_file
		sed -i -e "/${protein}/d" $sabanci_query_file
		sed -i -e "/${protein}/d" $truba_query_file
                let "count_completed=count_completed+1"
        elif grep -Fxq $protein $sabanci_query_file  ||  grep -Fxq $protein $truba_query_file ; then
                let "count_under_computation=count_under_computation+1"
        else
        	echo "$protein" >> $not_started_yet_query_file
        	let "count_incompleted=count_incompleted+1"
        fi
done
num_all=`cat $all_query_file|wc -l`
num_under_sabanci=`cat $sabanci_query_file|wc -l`
num_under_truba=`cat $truba_query_file|wc -l`

echo 
echo "All protein num: $num_all"
echo "Completed protein num: $count_completed"
echo "Not computed yet protein num: $count_incompleted"
echo "Ongoing computation protein num under truba: $num_under_truba, total: $count_under_computation"
echo "Ongoing computation protein num under truba: $num_under_sabanci, total: $count_under_computation"
