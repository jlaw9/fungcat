#source virtual-envs/py3env/bin/activate
python="python"
# try using numba for dot product which only works in anaconda
#python="/data/jeff-law/tools/anaconda3/bin/python"

#walltime="30:00:00"
declare -a version_list=("2018_09-s200-seq-sim-e0_1")
#cores="4"
cores="3"
#cores="12"
#walltime="120:00:00"
walltime="240:00:00"
declare -a algorithms=(
#"-A birgrank --alpha 0.5 --theta 0.5 --mu 0.5 --eps 0.0001"
"-A birgrank --alpha 0.75 --theta 0.5 --mu 0.5 --br-lambda 0.01 --eps 0.0001"
#"-A birgrank --alpha 0.9 --theta 0.5 --mu 0.5 --br-lambda 0.01 --eps 0.0001"
#"-A aptrank --diff-type twoway"
#"-A aptrank --diff-type oneway"
)
declare -a alg_names=("birgrank")
#declare -a alg_names=("aptrank-twoway" "aptrank-oneway")
#declare -a alg_names=("aptrank-twoway")

# I used squeeze as the dir for the alpha and num iterations evaluations, and eval-weights for everything else
#comparison="eval-weights"
comparison="eval-loso"
#comparison="squeeze"


## Evaluate cross-validation or LOSO of experimental annotations
#ev_codes="expc-rem-neg-comp-iea"
# uncomment these two for running the 200 species
#ann_cutoff="50" 
date="2018_09/"

# Evaluation 3: predict EXP+COMP, eval COMP
#ev_codes="expc-comp-rem-neg-iea"
ev_codes="expc-comp"
ann_cutoff="50" 
#eval_ev_codes="comp"  # a
#eval_ev_codes_cutoff="20"

# Evaluation 2: predict EXP+COMP, eval IEA
eval_ev_codes="iea"
eval_ev_codes_cutoff="50"

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
# for birgrank, parallelize the species because all GO term predictions are made simultaneously for each protein
for ((curr_idx = 1; curr_idx <= 40; curr_idx ++)); do
# finish the ones that stopped early
#declare -a iter_list=("11" "12")
#for curr_idx in ${iter_list[@]}; do
    alg=${algorithms[$i]}
    alg_name=${alg_names[$i]}

    exp_name="${th}${ev_codes}-50-1000${eval_ev_codes}${string}-${h}"
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
    # idx for selecting the set of species to run
    strains_file="inputs/selected-strains/s200-rand-groups/group-$curr_idx.txt"
    log_file="log/$comparison/${th}${ev_codes}${eval_ev_codes}/parallel-${alg_name}/$version-$real_exp_name-${alg_name}${weight_str}-$curr_idx.txt"
    mkdir -p log/$comparison/${th}${ev_codes}${eval_ev_codes}/parallel-${alg_name}

    # if the only_functions_file doesn't exist, then don't run it
    if [ ! -f "$strains_file" ]; then
        echo -e "\t strains_file doesn't exist. skipping: $strains_file"
    else
        # get the taxons
        taxons="" 
        for taxon in `cut -f 1 $strains_file`; do taxons="$taxons -T $taxon"; done; 
        echo $taxons
        cmd="""time $python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
            --version $version \
            --exp-name ${exp_name} \
            --pos-neg-file inputs/pos-neg/${date}${th_date}$ev_codes/pos-neg-${h}-${ann_cutoff}-list.tsv \
            $pos_neg_file_eval \
            -W 0  $alg $use_neg_opt $keep_ann_opt $weight_option $string_option \
            $taxons \
            --num-test-cutoff 10 \
            --postfix -group${curr_idx} \
            >> $log_file 2>&1
        """
            #--forcealg \
        echo $cmd
        # if test is specified, then don't actually run anything
        if [ "$1" == "--test" ]; then exit; fi
        echo $cmd >> $log_file

        # setup the qsub file
        echo "#PBS -l nodes=1:ppn=$cores,walltime=$walltime" > $qsub_file
        echo "#PBS -N $version-$exp_name-${alg_names[$i]}" >> $qsub_file
        # set the output files
        echo "#PBS -o $log_file-qsub.log" >> $qsub_file
        echo "#PBS -e $log_file-qsub.log" >> $qsub_file
        echo "" >> $qsub_file
        echo "cd $PWD" >> $qsub_file
        echo "echo \"Job started at \`date\`\" >> $log_file " >> $qsub_file
        if [ "$python" == "python" ]; then 
            echo "source virtual-envs/py3env-baobab/bin/activate" >> $qsub_file
        fi
        echo $cmd >> $qsub_file
        echo "echo \"Job finished at \`date\`\" >> $log_file " >> $qsub_file

        qsub $qsub_file
        sleep 10
    fi
done
done
done
done


