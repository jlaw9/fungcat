# setup for birgrank
# use the regular python for timing the algorithms to make sure there's no extra parallelism
#python="python"
# use anaconda so the scripts can be run on any ubuntu machine
python="/data/jeff-law/tools/anaconda3/bin/python"

version="2018_06-seq-sim-e0_1"
# for MF, don't use STRING
declare -a h_list=("mf")
string="-core"
if [ "$string" != "" ]; then
    version="${version}-string-700"
    string_option="--weight-swsn --string$string"
    weight_str="-swsn"
    declare -a h_list=("bp")
fi
cores="2"
walltime="20:00:00"
#ev_codes="expc-rem-neg-comp-iea"
ev_codes="expc"
only_func_ev_codes="expc"
ann_cutoff="10" 
#exp_name="eval-species-${ev_codes}-50-1000-wo-bug"
use_neg="-use-neg"
# this is the dir for the log files
comparison="eval-birgrank"

#declare -a alpha_list=("0.5")
declare -a alpha_list=("0.99" "0.95" "0.9" "0.75" "0.25" "0.1" "0.01")
#declare -a alpha_list=("0.9")
#declare -a theta_list=("1" "0.5")
declare -a theta_list=("0.5")
declare -a mu_list=("0.5")
#declare -a mu_list=("0.01" "0.99" "0.1" "0.9" "0.25" "0.75")
declare -a lambda_list=("0.5")
#declare -a lambda_list=("0.99" "0.9" "0.75" "0.25" "0.1" "0.01")
#declare -a lambda_list=("0.01")
declare -a eps_list=("0.0001")
#declare -a eps_list=("1e-5" "1e-6")
for alpha in ${alpha_list[@]}; do
for theta in ${theta_list[@]}; do
for mu in ${mu_list[@]}; do
for lambda in ${lambda_list[@]}; do
for eps in ${eps_list[@]}; do
for h in ${h_list[@]}; do

    exp_name="${ev_codes}-50-1000${eval_ev_codes}${string}-$h"
    real_exp_name="${exp_name}$use_neg$keep_ann"
    log_dir="log/$comparison/$real_exp_name/$version"
    log_file="$log_dir/birgrank-eps$eps-a$alpha-t$theta-m$mu-l$lambda${weight_str}.txt"
    mkdir -p $log_dir
    qsub_file="$PWD/$log_dir/birgrank-eps$eps-a$alpha-t$theta-m$mu-l$lambda${weight_str}.qsub"

    if [ "$eval_ev_codes" != "" ]; then
    #if [ "$th_date" != "" ]; then
        pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/${date}/$eval_ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
    fi
    only_functions_file="inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000.txt"

    cmd="""time $python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
		--version $version \
		--exp-name ${exp_name} \
		--pos-neg-file inputs/pos-neg/${date}${th_date}$ev_codes/pos-neg-${h}-${ann_cutoff}-list.tsv \
		--only-functions $only_functions_file \
        $pos_neg_file_eval \
		-W 0  -A birgrank \
        --alpha $alpha --theta $theta --mu $mu --br-lambda $lambda $use_neg_opt $keep_ann_opt $string_option \
        --eps $eps \
        --num-test-cutoff 10 \
        >> $log_file 2>&1
    """
		#--forcealg \
    echo $cmd
    # if test is specified, then don't actually run anything
    if [ "$1" == "--test" ]; then
        exit
    fi
    echo $cmd >> $log_file

    # setup the qsub file
    echo "#PBS -l nodes=1:ppn=$cores,walltime=$walltime" > $qsub_file
    echo "#PBS -N $version-$exp_name-a$alpha-t$theta-m$mu" >> $qsub_file
    # set the output files
    echo "#PBS -o $qsub_file-out.log" >> $qsub_file
    echo "#PBS -e $qsub_file-out.log" >> $qsub_file
    echo "" >> $qsub_file
    echo "cd $PWD" >> $qsub_file
    echo "echo \"Job started at \`date\`\" >> $log_file " >> $qsub_file
    if [ "$python" == "python" ]; then
        echo "source virtual-envs/py3env-baobab/bin/activate" >> $qsub_file
    fi
    echo $cmd >> $qsub_file
    echo "echo \"Job finished at \`date\`\" >> $log_file " >> $qsub_file

    #machine="baobab"
    #machine="spode"
    #machine="agatha"
    #machine="cuthbert"
    #machine="prosser"
    #machine="wyatt"
    #machine="molloy"
    #machine="honoria"
    machine="pirbright"
    #machine="cowcreamer"
    #machine="simmons"
    # can't start jobs on these. screen version is too old
    #machine="carmody"
    #machine="mnemosyne"

    screen_name="$exp_name-birgrank-eps$eps-a$alpha-t$theta-m$mu-l$lambda-$version" 
    if [ "$machine" == "baobab" ]; then
        if [ "`hostname`" == "baobab.cbb.lan" ]; then
            qsub $qsub_file
        elif [ "`hostname`" == "$machine" ]; then
            ssh $machine "qsub $qsub_file"
        fi
    # if this is the machine we're currently on, then no need to ssh
    elif [ "`hostname`" == "$machine" ]; then
        echo """screen -S $screen_name -d -m /bin/sh -c \"bash $qsub_file\""""
        screen -S $screen_name -d -m /bin/sh -c "bash $qsub_file"
    else
        # the qsub file has the full file path, so we can just use that
        echo """
        ssh $machine \"screen -S $screen_name -d -m /bin/sh -c \
            \\\"bash $qsub_file\\\"\"
        """
        ssh $machine "screen -S $screen_name -d -m /bin/sh -c \
            \"bash $qsub_file\""
    fi

done
done
done 
done
done
done

