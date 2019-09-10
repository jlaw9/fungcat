source virtual-envs/py3env-ubuntu14/bin/activate

version="2018_06-seq-sim-e0_1"
#version="2018_09-s200-seq-sim-e0_1"
#version="2018_06-seq-sim-e0_1-string"
ev_codes="expc-rem-neg-comp-iea"
exp_name="${ev_codes}"
only_functions="inputs/only-functions/expc/expc-50-1000.txt"
#exp_name="${ev_codes}-50-1000"
alg="sinksource"
#weight="--weight-per-goterm --string-core"

#declare -a alpha_list=("1" "0.9" "0.8" "0.7" "0.6" "0.5")
#declare -a alpha_list=("0.99" "0.95" "0.9" "0.8" "0.7" "0.6" "0.5")
declare -a alpha_list=("1" "0.99" "0.975")
max_iters=1000
comparison='alpha'
#declare -a max_iters_list=("400" "200" "50" "20" "10")
#declare -a max_iters_list=("5" "2" "1")
#declare -a max_iters_list=("20")
#alpha=0.8
#comparison='max-iters'
for alpha in ${alpha_list[@]}; do
#for max_iters in ${max_iters_list[@]}; do
#    exp_name="${ev_codes}-50-1000-iters${max_iters}"

    #max_iters="" 
    #if [ "$alpha" == "1.0" ]; then
    #    max_iters="--max-iters 1000 "
    #fi
    log_file="log/squeeze/compare-$comparison/${version}-${exp_name}-${alg}-${alpha}.log"

    cmd="""
    time python -u src/algorithms/gain_scipy/run_algs.py \
        --version $version \
        --exp-name $exp_name \
        --pos-neg-file inputs/pos-neg/$ev_codes/pos-neg-mf-10-list.tsv \
        --only-functions $only_functions  \
        -A $alg $weight \
        -a $alpha \
        --forcealg \
        --num-pred-to-write 0 \
        --eps 0  --max-iters $max_iters \
        --only-cv -C 5 \
        >> $log_file 2>&1
    """
        #--pos-neg-file inputs/pos-neg/$ev_codes/pos-neg-bp-10-list.tsv \
        #--eps 0.0001 \
        #-A sinksource-squeeze \
        #--epsUB 0.01 \
    echo "$cmd"
    echo "$cmd" >> $log_file
    # apparently you need to redirect the output here as well
    $cmd >> $log_file 2>&1

done
