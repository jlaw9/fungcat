#source virtual-envs/py3env/bin/activate
#python="python"
python="/data/jeff-law/tools/anaconda3/bin/python"

machine="baobab"
#machine="spode"
#machine="agatha"
#machine="cuthbert"
#machine="prosser"
#machine="parsloe"
#machine="wyatt"
#machine="molloy"
#machine="honoria"
#machine="pirbright"
#machine="cowcreamer"
#machine="simmons"
# can't start jobs on these. screen version is too old
#machine="carmody"
#machine="mnemosyne"

declare -a version_list=(
"2018_06-seq-sim-e0_1"
#"2018_06-seq-sim-e0_1-string-700" 
# try different STRING cutoffs 
#"2018_06-seq-sim-e0_1-string"
#"2018_06-seq-sim-e0_1-string-150" "2018_06-seq-sim-e0_1-string-900"
# try different BLAST E-value cutoffs
#"2018_06-seq-sim-e1e-6" "2018_06-seq-sim-e1e-15" "2018_06-seq-sim-e1e-25"
#"2018_06-seq-sim-e5" "2018_06-seq-sim-e20" "2018_06-seq-sim-e50"
#"2018_09-s200-seq-sim-e0_1"
)
# If the version is a STRING version, uncomment this line
#string="-core"
declare -a h_list=(
#"bp" 
"mf"
)

#cores="2"
# use 24 cores for aptrank
cores="24"
walltime="30:00:00"
# even though genemania doesn't use alpha, the output file is still written with it...
declare -a algorithms=(
#"-A sinksource --alpha 1.0 --eps 0.0001" 
#"-A sinksource --max-iters 20 --alpha 0.95 --eps 0.0" 
"-A sinksource --max-iters 10 --alpha 0.95 --eps 0.0" 
#"-A local -A sinksourceplus --max-iters 20 --alpha 0.95 --eps 0.0" 
#"-A genemania -A localplus --max-iters 20 --alpha 0.95 --eps 0.0"
#"-A birgrank --alpha 0.9 --theta 0.5 --mu 0.5 --br-lambda 0.01 --eps 0.0001"
#"-A localplus --max-iters 20 --alpha 0.95 --eps 0.0"
#"-A aptrank --alpha 0.9 --diff-type oneway --br-lambda 0.01"
#"-A sinksource-squeeze --alpha 0.95 --max-iters 1000 --eps 0.0 --rank-pos-neg --compare-ranks" 
)
declare -a alg_names=(
#"sinksource-a0.95-maxi20" 
"sinksource-a0.95-maxi10" 
#"sinksourceplus" 
#"genemania" 
#"birgrank"
#"localplus"
#"aptrank-oneway"
#"sinksource-squeeze-a0.95"
)

# I used eval-cutoffs for everything except speed-test
comparison="eval-loso"
#comparison="speed-test"
# use a single core to ensure its a fair comparison across the board
if [ "$comparison" == "speed-test" ]; then cores="1"; fi

# TODO construct this with a for loop or something
#declare -a algorithms=(
#"-A sinksource --max-iters 1000 --alpha 1.0 --eps 0.0" 
#"-A sinksource --max-iters 1000 --alpha 0.99 --eps 0.0"   
###"-A sinksource --max-iters 1000 --alpha 0.975 --eps 0.0"
#"-A sinksource --max-iters 1000 --alpha 0.95 --eps 0.0"  
#"-A sinksource --max-iters 1000 --alpha 0.9 --eps 0.0"
#"-A sinksource --max-iters 1000 --alpha 0.8 --eps 0.0" 
#"-A sinksource --max-iters 1000 --alpha 0.7 --eps 0.0" 
#"-A sinksource --max-iters 1000 --alpha 0.6 --eps 0.0" 
#"-A sinksource --max-iters 1000 --alpha 0.5 --eps 0.0" 
#"-A sinksource --max-iters 400 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 200 --alpha 0.95 --eps 0.0"
#"-A sinksource --max-iters 50 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 20 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 10 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 5 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 2 --alpha 0.95 --eps 0.0" 
#"-A sinksource --max-iters 1 --alpha 0.95 --eps 0.0" 
## don't compare ranks to time the algorithm
##"-A sinksource-squeeze --alpha 0.99 --max-iters 5000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.99 --max-iters 1000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.8 --max-iters 1000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.7 --max-iters 1000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.6 --max-iters 1000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.5 --max-iters 1000 --eps 0.0 --rank-pos-neg" 
##"-A sinksource-squeeze --alpha 0.95 --max-iters 1000 --eps 0.0 --rank-pos-neg --compare-ranks" 
##"-A sinksource-squeeze --alpha 0.9 --max-iters 1000 --eps 0.0 --rank-pos-neg --compare-ranks" 
##"-A genemania --max-iters 20 --alpha 0.95 --eps 0.0 --tol 1e-4"
#)
#declare -a alg_names=(
#"sinksource-a1.0-maxi1000"
#"sinksource-a0.99-maxi1000"  
##"sinksource-a0.975-maxi1000"  # also test these for MF
#"sinksource-a0.95-maxi1000" 
#"sinksource-a0.9-maxi1000"
#"sinksource-a0.8-maxi1000"
#"sinksource-a0.7-maxi1000"
#"sinksource-a0.6-maxi1000"
#"sinksource-a0.5-maxi1000"
#"sinksource-a0.95-maxi400"
#"sinksource-a0.95-maxi200"
#"sinksource-a0.95-maxi50"
#"sinksource-a0.95-maxi20"
#"sinksource-a0.95-maxi10"
#"sinksource-a0.95-maxi5"
#"sinksource-a0.95-maxi2"
#"sinksource-a0.95-maxi1"
##"sinksource-squeeze-a0.99-maxi5000"
##"sinksource-squeeze-a0.99"
##"sinksource-squeeze-a0.8"
##"sinksource-squeeze-a0.7"
##"sinksource-squeeze-a0.6"
##"sinksource-squeeze-a0.5"
##"sinksource-squeeze-a0.95"
##"sinksource-squeeze-a0.90"
##"genemania-tol1e-4"
#)

## Weight and STRING options
#weight_option="--unweighted"
#weight_str="-unw"
## option to unweight the STRING networks, but still combine the weights using the weight-per-goterm method
#weight_option="--unweighted --weight-gm2008"
#weight_str="-unw-gm2008"
# part of the experiment name
# currently the string networks are hard coded
#string="-non-transferred"
#string="-all"
if [ "$string" != "" ]; then 
    string_option="--string$string"
    ## option to integrate network weights using the geneMANIA findKernelWeights method per goterm
    #weight_option="--weight-gm2008"
    #weight_str="-gm2008"
    ## option to weight and combine the string and seq-seim networks once using all GO term annotations of a specific hierarchy. 
    ## needed for birgrank with STRING
    weight_option="--weight-swsn"
    weight_str="-swsn"
fi


## Evaluate cross-validation or LOSO of experimental annotations
#ev_codes="expc"
##ev_codes="expc-rem-neg-comp-iea"
#only_func_ev_codes="expc"
#ann_cutoff="10" 

# Evaluate EXP+COMP LOSO
ev_codes="expc-comp"
#ev_codes="expc-comp-rem-neg-iea"
only_func_ev_codes="expc-comp"
ann_cutoff="10" 

# Evaluation 3: predict EXP+COMP, eval COMP or COMP+IEA
#eval_ev_codes="comp"  # a
#eval_ev_codes_cutoff="5"
##eval_ev_codes="comp-iea"  # b
eval_ev_codes="iea"  # c
eval_ev_codes_cutoff="10"

# for 200 species
# Evaluation: predict EXP+COMP, eval COMP
#ann_cutoff="50" 
#eval_ev_codes="comp"  # a
#eval_ev_codes_cutoff="20"
# Evaluation 2: predict EXP+COMP, eval IEA
#eval_ev_codes="iea"
#only_func_ev_codes="iea"
#eval_ev_codes_cutoff="50"
# UPDATE: for recover IEA, switch off between predicting all species at once for species which have no left-out annotations,
# and leaving out species which have annotations to be left-out
#if [ "$eval_ev_codes" == "iea" ]; then
#    #keep_ann="-keep-ann"
#    #keep_ann_opt="--keep-ann"
#    keep_ann=""
#    keep_ann_opt="--eval-goterms-with-left-out-only"
#fi


# Temporal Holdout evaluation
#th="th-"
#th_date="2016_06/"
#ev_codes="expc-rem-neg-comp-iea"
##eval_ev_codes="expc-rem-neg-comp-iea"
#eval_ev_codes_cutoff="10"
#only_func_ev_codes="expc"
#ann_cutoff="10" 
##keep_ann="-keep-ann"
##keep_ann_opt="--keep-ann"
#keep_ann="-oracle-keep-ann"
#keep_ann_opt="--oracle --keep-ann"

# also compare with/without using negative examples for evaluation
# only used for the leave-one-species-out evaluation
use_neg="-use-neg"
#use_neg="-non-pos-neg"
#use_neg_opt="--non-pos-as-neg-eval"

# for the rest of the algorithms
for h in ${h_list[@]}; do
for version in ${version_list[@]}; do
#for alg in ${algorithms[@]}; do
for ((i = 0; i < ${#algorithms[@]}; i++)); do
    alg=${algorithms[$i]}

    exp_name="${th}${ev_codes}-50-1000${eval_ev_codes}${string}-${h}"
    if [ "$comparison" == "speed-test" ]; then
        exp_name="speed-$exp_name"
    fi
    # it's really 50+ for s200, but it's too much of a hassle to change in other scripts
    #if [ "$version" == "2018_09-s200-seq-sim-e0_1" ]; then
    #    exp_name="${th}${ev_codes}-50${eval_ev_codes}${string}-${h}"
    #fi
    #exp_name="parallel-test-${h}"
    real_exp_name="${exp_name}$use_neg$keep_ann"
    log_dir="log/$comparison/$real_exp_name/$version"
    log_file="$log_dir/${alg_names[$i]}${weight_str}.txt"
    mkdir -p $log_dir
    qsub_file="$PWD/$log_dir/${alg_names[$i]}${weight_str}.qsub"

    if [ "$eval_ev_codes" != "" ]; then
    #if [ "$th_date" != "" ]; then
        pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/${date}/$eval_ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
        #pos_neg_file_eval="--pos-neg-file-eval inputs/pos-neg/$ev_codes/pos-neg-${h}-${eval_ev_codes_cutoff}-list.tsv"
    fi
    #--only-functions inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000.txt \
    only_functions_file="inputs/only-functions/$only_func_ev_codes/${only_func_ev_codes}-50-1000.txt"

    # if the only_functions_file doesn't exist, then don't run it
    if [ ! -f "$only_functions_file" ]; then
        echo -e "\t only_functions_file doesn't exist. skipping: $only_functions_file"
    else
        # swap out the evaluation to run
        cmd="""time $python -u src/algorithms/gain_scipy/run_algs.py --only-cv -C 5 \
            """
        cmd="""time $python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
            --version $version \
            --exp-name ${exp_name} \
            --pos-neg-file inputs/pos-neg/${date}${th_date}$ev_codes/pos-neg-${h}-${ann_cutoff}-list.tsv \
            --only-functions $only_functions_file \
            $pos_neg_file_eval \
            -W 0  $alg $use_neg_opt $keep_ann_opt $weight_option $string_option \
            --num-test-cutoff 10 \
            --forcealg \
            >> $log_file 2>&1
        """
        echo $cmd
        # if test is specified, then don't actually run anything
        if [ "$1" == "--test" ]; then exit; fi
        echo $cmd >> $log_file

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

        #qsub $qsub_file
        # this will run the command instead of submit it to the queue
#        $cmd >> $log_file

        screen_name="$version-$exp_name-${alg_names[$i]}" 
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
    fi
done
done
done

