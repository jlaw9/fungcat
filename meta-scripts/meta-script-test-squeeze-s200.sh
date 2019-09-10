#source virtual-envs/py3env/bin/activate
# anaconda has some parallelization that my regular python environment doesn't have. When I'm timing the algorithms, I want to run it without parallelization
python="/data/jeff-law/tools/anaconda3/bin/python"

# so far I've been running things here rather than on baobab
# so I can make sure each computer only has one thing running on it
# for a more fair time comparison

machine="cowcreamer"
#machine="honoria"
#machine="parsloe"
#machine="pirbright"

version="2018_09-s200-seq-sim-e0_1"
date="2018_09/"
#ev_codes="expc-comp-rem-neg-iea"
ev_codes="expc-comp"
#eval_ev_codes="comp"
#eval_ev_codes_cutoff="20"
#only_functions="inputs/only-functions/${date}/comp/comp-${eval_ev_codes_cutoff}-rand-50.txt"

eval_ev_codes="iea"
eval_ev_codes_cutoff="50"
only_functions="inputs/only-functions/${date}/iea/iea-${eval_ev_codes_cutoff}-rand-50.txt"
# for recover IEA, switch off between predicting all species at once for species which have no left-out annotations,
# and leaving out species which have annotations to be left-out
if [ "$eval_ev_codes" == "iea" ]; then
    #keep_ann="-keep-ann"
    #keep_ann_opt="--keep-ann"
    keep_ann=""
    keep_ann_opt="--eval-goterms-with-left-out-only"
fi

alg="sinksource-squeeze"
comparison='compare-cutoffs'
#comparison='ranks-all'
#options="--rank-all"
options="--rank-pos-neg --compare-ranks"
#alg="sinksource"
#options="--only-cv -C 5"
#declare -a alpha_list=("1" "0.9" "0.8" "0.7" "0.6" "0.5")
#declare -a alpha_list=("0.9" "0.8" "0.7" "0.6" "0.5")
# for MF, I tested more alpha values
#declare -a alpha_list=("1" "0.99" "0.95" "0.9" "0.8" "0.7" "0.6" "0.5")
#declare -a alpha_list=("0.9" "0.95" "0.99" "0.8" "0.7" "0.6" "0.5")
#declare -a alpha_list=("1")
#declare -a alpha_list=("0.99")
declare -a alpha_list=("0.95")
#declare -a alpha_list=("0.8")
#comparison='alpha'
#declare -a max_iters_list=("50" "20" "10")
#declare -a max_iters_list=("400" "200" "50" "20" "10" "5" "2" "1")
#declare -a max_iters_list=("1" "2" "5" "10" "20" "50" "200" "400")
#declare -a max_iters_list=("5" "2" "1")
#declare -a max_iters_list=("20")
declare -a max_iters_list=("5000")
#declare -a max_iters_list=("0")  # run conjugate gradient
#comparison='max-iters'
#declare -a h_list=("bp" "mf")
declare -a h_list=("bp")
#declare -a h_list=("mf")
for h in ${h_list[@]}; do
for alpha in ${alpha_list[@]}; do
for max_iters in ${max_iters_list[@]}; do
    exp_name="${ev_codes}-${eval_ev_codes_cutoff}-rand-50${eval_ev_codes}${string}-${h}"
    real_exp_name="${exp_name}$use_neg$keep_ann"

    log_dir="log/$comparison/$real_exp_name/$version"
    log_file="$log_dir/$alg-a${alpha}-maxi${max_iters}.txt"
    mkdir -p $log_dir

    if [ "$eval_ev_codes" != "" ]; then
    #if [ "$th_date" != "" ]; then
        pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/${date}/$eval_ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
    fi
    cmd="""time $python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
        --version $version  \
        --exp-name $exp_name \
        --pos-neg-file inputs/pos-neg/${date}/$ev_codes/pos-neg-${h}-50-list.tsv \
        $pos_neg_file_eval \
        --only-functions $only_functions  \
        -A $alg $options \
        -a $alpha \
        --forcealg \
        --num-pred-to-write 0 --num-test-cutoff 10 $keep_ann_opt \
        --eps 0  --max-iters $max_iters \
        >> $log_file 2>&1
    """
        #--eps 0.0001 \
        #--epsUB 0.01 \
    echo "$cmd"
    # if test is specified, then don't actually run anything
    if [ "$1" == "--test" ]; then exit; fi
    echo "$cmd" >> $log_file
    # apparently you need to redirect the output here as well
    #$cmd >> $log_file 2>&1
    cmd="cd $PWD; $cmd"

    screen_name="$version-$exp_name-$alg" 
    #if [ "$machine" == "baobab" ]; then
    #    if [ "`hostname`" == "baobab.cbb.lan" ]; then
    #        qsub $qsub_file
    #    elif [ "`hostname`" == "$machine" ]; then
    #        ssh $machine "qsub $qsub_file"
    #    fi
    # if this is the machine we're currently on, then no need to ssh
    if [ "`hostname`" == "$machine" ]; then
        echo """screen -S $screen_name -d -m /bin/sh -c \"$cmd\""""
        screen -S $screen_name -d -m /bin/sh -c "$cmd"
    else
        echo """
        ssh $machine \"screen -S $screen_name -d -m /bin/sh -c \
            \\\"$cmd\\\"\"
        """
        ssh $machine "screen -S $screen_name -d -m /bin/sh -c \
            \"$cmd\""
    fi

done
done
done


