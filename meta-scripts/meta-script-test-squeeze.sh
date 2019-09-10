#source virtual-envs/py3env-ubuntu14/bin/activate

# so far I've been running things here rather than on baobab
# so I can make sure each computer only has one thing running on it
# for a more fair time comparison

string="-core"
version="2018_06-seq-sim-e0_1"
#version="2018_09-s200-seq-sim-e0_1"
if [ "$string" != "" ]; then
    version="${version}-string-700"
    #weight="--weight-gm2008 --string-core"
    weight="--weight-swsn --string-core"
fi
ev_codes="expc"
#ev_codes="expc-rem-neg-comp-iea"
#exp_name="${ev_codes}"
#only_functions="inputs/only-functions/expc/expc-50-1000.txt"
only_functions="inputs/only-functions/expc/expc-50-1000-rand-0_05/expc-50-1000-rand-0_05.txt"
#only_functions="inputs/only-functions/2018_09/expc/expc-50-rand-50.txt"

alg="sinksource-squeeze"
comparison='ranks'
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
declare -a max_iters_list=("1000")
#declare -a max_iters_list=("0")  # run conjugate gradient
#comparison='max-iters'
#declare -a h_list=("bp" "mf")
declare -a h_list=("bp")
#declare -a h_list=("mf")
for h in ${h_list[@]}; do
for alpha in ${alpha_list[@]}; do
for max_iters in ${max_iters_list[@]}; do
    #exp_name="${ev_codes}-50-1000${string}-${h}"
    exp_name="${ev_codes}-50-1000${string}-${h}-rand0_05"
      

    #max_iters="" 
    #if [ "$alpha" == "1.0" ]; then
    #    max_iters="--max-iters 1000 "
    #fi
    log_file="log/squeeze/compare-$comparison/${version}-${exp_name}-${alg}-a${alpha}-maxi${max_iters}.log"
    mkdir -p log/squeeze/compare-$comparison

    #cmd="""time python -u src/algorithms/gain_scipy/run_algs.py \
    #    """
    # anaconda has some parallelization that my regular python environment doesn't have. When I'm timing the algorithms, I want to run it without parallelization
    #cmd="""time /data/jeff-law/tools/anaconda3/bin/python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
    cmd="""time python -u src/algorithms/gain_scipy/eval_leave_one_species_out.py \
        --version $version $weight \
        --exp-name $exp_name \
        --pos-neg-file inputs/pos-neg/$ev_codes/pos-neg-${h}-10-list.tsv \
        --only-functions $only_functions  \
        -A $alg $options \
        -a $alpha \
        --forcealg \
        --num-pred-to-write 0 --num-test-cutoff 10 \
        --eps 0  --max-iters $max_iters \
        >> $log_file 2>&1
    """
        #--only-cv -C 5 \
        #--eps 0.0001 \
        #--epsUB 0.01 \
    echo "$cmd"
    echo "$cmd" >> $log_file
    # apparently you need to redirect the output here as well
    $cmd >> $log_file 2>&1

done
done
done

