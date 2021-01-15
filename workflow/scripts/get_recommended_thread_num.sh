cat $1|grep MPI|awk -F: '{print $2}'
