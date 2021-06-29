#./select_n_query.sh sabanci 512
cluster_name=$1 # sabanci or truba
num_query=$2

not_started_yet_query_file=../not_yet_computed
sabanci_query_file=../under_computation_tosun
truba_query_file=../under_computation_levrek1

num_remain_query=`cat $not_started_yet_query_file|wc -l`
if [ "$num_query" -gt "$num_remain_query" ]; then
	echo "max. $num_remain_query query set to run"
	exit 1
fi

if [ "$cluster_name" == "sabanci" ]; then
	head -n $num_query $not_started_yet_query_file > $sabanci_query_file
	sed -i -e "1,${num_query}d" $not_started_yet_query_file
elif [ "$cluster_name" == "truba" ]; then
	head -n $num_query $not_started_yet_query_file > $truba_query_file
	sed -i -e "1,${num_query}d" $not_started_yet_query_file
else
	echo "select truba or sabanci cluster"
fi
