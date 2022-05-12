###############################################################
# Trap signals
###############################################################
notify_term () {
    echo "... remove jobname $1"
    scancel --name $1
    exit 9
}

###############################################################
# Lock functions
###############################################################
pi_lock () {
    if [ "x1" == "x"$(find . -name "${1//\.[A-Z0-9]*/}*.lock" -printf 1) ]; then
        echo >&2 "ERROR process with jobname $1 can't be started as a lock is present"
        kill -s INT $$
    elif [ -e ${1//\.[A-Z0-9]*/}.done ]; then
        echo 0 
    else
        echo 0 > ${1}.lock
        echo 1
    fi
}

pi_unlock () {
    if [ -e ${1}.lock ]; then
        if [ "x0" == "x"$(jobs_failures ${1}) ]; then
            mv ${1}.lock ${1//\.[A-Z0-9]*/}.done
            echo_info "Unlock jobname $1"       
        else
            echo >&2 "ERROR in execution of jobname $1"
            exit 1
        fi
    fi
}

pi_wait () {
    trap "notify_term $1" 1 2 3 8 9 14 15
    sleep 10 # remember there are latencies with filesystems
    if [ -e ${1}.lock ]; then
        # ensure the job started before the real waiting
        if [ "x0" == "x"$(head -1 ${1}.lock) ]; then
            if [ "x0" == "x"$(jobs_total ${1}) ]; then
                pi_wait ${1}
            fi 
            echo $$ > ${1}.lock
        fi
            
        # wait pending jobs
        n=$(jobs_pending ${1})
        if [ "x0" == "x$n" ]; then
            pi_unlock ${1}
        else
            echo "... waiting $n jobs for task ${1}"
            pi_wait ${1}
        fi
    fi
}

jobs_pending () {
    n=$(sacct --name ${1} -n -b -p | egrep "RUNNING|PENDING" | wc -l)
    echo $n
}

jobs_failures () {
    n=$(sacct --name ${1} -n -b -p  | grep -v 'COMPLETED|0:0' | wc -l)
    echo $n
}

jobs_total () {
    n=$(sacct --name ${1} -n -b -p | wc -l)
    echo $n
}

###############################################################
# UID generator
###############################################################
uid_gen () {
    < /dev/urandom tr -dc A-Z0-9 | head -c${1:-8}
    echo
}

###############################################################
# Prepare dir
###############################################################
set_dir () {
    if [ "x" != "x${1}" ]; then
        if [ -d ${1} ]; then
            rm -rf ${1};
            rm -rf ${1}_logs
        fi
        mkdir ${1}
        mkdir ${1}_logs
    fi
}

###############################################################
# Show info with date
###############################################################
echo_info () {
    d=$(date)
    echo "[$d] $1"
}

###############################################################
# Force value to string
###############################################################
to_string () {
    printf '%s' $1    
}

