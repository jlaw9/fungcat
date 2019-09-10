#source virtual-envs/py3env/bin/activate
python="python"
# try using numba for dot product which only works in anaconda
#python="/data/jeff-law/tools/anaconda3/bin/python"

cores="1"
#cores="3"
declare -a version_list=("2018_09-s200-seq-sim-e0_1")
#cores="6"
#cores="12"
#walltime="120:00:00"
walltime="240:00:00"
# even though genemania doesn't use alpha, the output file is still written with it...
declare -a algorithms=(
#"-A genemania --max-iters 20 --alpha 0.95 --eps 0.0"
#"-A birgrank --alpha 0.9 --theta 0.5 --mu 0.5 --br-lambda 0.01 --eps 0.0001"
"-A sinksource --max-iters 10 --alpha 0.95 --eps 0.0" 
#"-A localplus --max-iters 20 --alpha 0.95 --eps 0.0" 
#"-A local -A sinksourceplus --max-iters 20 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 1000 --alpha 1 --eps 0.0" 
#"-A sinksource --max-iters 50 --alpha 1 --eps 0.0" 
#"-A sinksource --max-iters 10 --alpha 1 --eps 0.0" 
#"-A sinksource --max-iters 20 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 50 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 200 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 400 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 5 --alpha 0.95 --eps 0.0" 
)
declare -a alg_names=(
#"genemania"
#"birgrank"
"sinksource-a0.95-maxi10"
#"localplus"
#"sinksourceplus"
#"sinksource-a1-maxi1000"
#"sinksource-a1-maxi50"
#"sinksource-a1-maxi10"
#"sinksource-a0.95-maxi20"
#"sinksource-a0.95-maxi50"
#"sinksource-a0.95-maxi200"
#"sinksource-a0.95-maxi400"
#"sinksource-a0.95-maxi5"
)

# I used squeeze as the dir for the alpha and num iterations evaluations, and eval-loso for everything else
#comparison="eval-loso"
comparison="speed-test"
#comparison="squeeze"


## Evaluate cross-validation or LOSO of experimental annotations
#ev_codes="expc"
#only_func_ev_codes="expc"
# uncomment these two for running the 200 species
ann_cutoff="50" 
date="2018_09/"

# Evaluation 3: predict EXP+COMP, eval COMP
##ev_codes="expc-comp-rem-neg-iea"
ev_codes="expc-comp"
only_func_ev_codes="expc-comp"
ann_cutoff="50" 
eval_ev_codes="comp"  # a
eval_ev_codes_cutoff="20"

# Evaluation 2: predict EXP+COMP, eval IEA
#eval_ev_codes="iea"
#only_func_ev_codes="iea"
#eval_ev_codes_cutoff="50"
# UPDATE: for recover IEA, run both of these modes: predicting all species at once for those which have no left-out annotations (keep-ann),
# and leaving out species which have annotations to be left-out (eval goterms with left-out only)
if [ "$eval_ev_codes" == "iea" ]; then
    keep_ann="-keep-ann"
    keep_ann_opt="--keep-ann"
    #keep_ann=""
    #keep_ann_opt="--eval-goterms-with-left-out-only"
fi

# also compare with/without using negative examples for evaluation
# only used for the leave-one-species-out evaluation
use_neg="-use-neg"
#use_neg="-non-pos-neg"
#use_neg_opt="--non-pos-as-neg-eval"

# for the rest of the algorithms
#declare -a h_list=("bp" "mf")
declare -a h_list=("bp")
#declare -a h_list=("mf")
for h in ${h_list[@]}; do
for version in ${version_list[@]}; do
#for alg in ${algorithms[@]}; do
for ((i = 0; i < ${#algorithms[@]}; i++)); do
# UPDATE: parallelize the GO terms
# with the jobs split into 48 random groups, each group takes about 9% of RAM of the 64GB on a baobab node, meaning I can run about 10 in parallel. 
# To split it up evenly and give each job a little cushion room for RAM, I give each job 3 cores, which would be 8 running in parallel
#for ((curr_idx = 1; curr_idx <= 48; curr_idx ++)); do
# finish the ones that stopped early
#declare -a iter_list=("1" "2" "3" "6" "7" "8" "16" "17" "18" "19" "27" "28" "29" "30" "40" "41" "48")
#for curr_idx in ${iter_list[@]}; do
    alg=${algorithms[$i]}
    alg_name=${alg_names[$i]}

    exp_name="${th}${ev_codes}-50-1000${eval_ev_codes}${string}-${h}"
    if [ "$comparison" == "speed-test" ]; then
        exp_name="speed-${th}${ev_codes}-50-1000${eval_ev_codes}${string}-${h}"
    fi
    # it's really 50+ for s200, but it's too much of a hassle to change in other scripts
    #if [ "$version" == "2018_09-s200-seq-sim-e0_1" ]; then
    #    exp_name="${th}${ev_codes}-50${eval_ev_codes}${string}-${h}"
    #fi
    #exp_name="parallel-test-${h}"
    real_exp_name="${exp_name}$use_neg$keep_ann"
    log_file="log/$comparison/${th}${ev_codes}${eval_ev_codes}/$version-$real_exp_name-${alg_names[$i]}${weight_str}.txt"
    mkdir -p log/$comparison/${th}${ev_codes}${eval_ev_codes}
    qsub_file="$PWD/log/$comparison/${th}${ev_codes}${eval_ev_codes}/$version-$real_exp_name-${alg_names[$i]}.qsub"

    if [ "$eval_ev_codes" != "" ]; then
    #if [ "$th_date" != "" ]; then
        pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/${date}/$eval_ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
        #pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/$ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
    fi
    # idx for selecting the only functions
    # TODO set the normal only_functions_file
    if [ "$curr_idx" != "" ]; then
        if [ "$version" == "2018_09-s200-seq-sim-e0_1" ]; then
            only_functions_file="inputs/only-functions/${date}/$only_func_ev_codes/${only_func_ev_codes}-rand-groups/group-$curr_idx.txt"
        else
            only_functions_file="inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000-rand-groups/group-$curr_idx.txt"
        fi
        log_file="log/$comparison/${th}${ev_codes}${eval_ev_codes}/parallel-${alg_name}/$version-$real_exp_name-${alg_names[$i]}${weight_str}-$curr_idx.txt"
        mkdir -p log/$comparison/${th}${ev_codes}${eval_ev_codes}/parallel-${alg_name}
    else
        #--only-functions inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000.txt \
        only_functions_file="inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000.txt"
    fi

    # if the only_functions_file doesn't exist, then don't run it
    if [ ! -f "$only_functions_file" ]; then
        echo -e "\t only_functions_file doesn't exist. skipping: $only_functions_file"
    #elif [ ! -f "$log_file" ]; then
    #    echo "skipping $log_file"
    else
        # swap out the evaluation to run
        cmd="""time $python -u src/algorithms/gain_scipy/run_algs.py --only-cv -C 5 \
            """
        cmd="""time $python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
            --version $version \
            --exp-name ${exp_name} \
            --pos-neg-file inputs/pos-neg/${date}${th_date}$ev_codes/pos-neg-${h}-${ann_cutoff}-list.tsv \
            $pos_neg_file_eval \
            -W 0  $alg $use_neg_opt $keep_ann_opt $weight_option $string_option \
            --only-functions $only_functions_file \
            --num-test-cutoff 10 \
            >> $log_file 2>&1
        """
            #--postfix -group${curr_idx} \
            #--forcealg \
        echo $cmd
        # if test is specified, then don't actually run anything
        if [ "$1" == "--test" ]; then exit; fi

        # setup the qsub file
        echo "#PBS -l nodes=1:ppn=$cores,walltime=$walltime" > $qsub_file
        echo "#PBS -N $version-$exp_name-${alg_names[$i]}" >> $qsub_file
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

        echo $cmd >> $log_file
        qsub $qsub_file
    fi
#done
done
done
done

