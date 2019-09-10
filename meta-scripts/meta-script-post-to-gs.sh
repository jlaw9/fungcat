source virtual-envs/py3env/bin/activate

#ev_codes="expc-rem-neg-comp-iea"
ev_codes="expc"
ann_cutoff="10"
#eval_ev_codes="comp"
eval_ev_codes_cutoff="5"
#exp_name="eval-sp-${ev_codes}-50-1000-${eval_ev_codes}-bp"
use_neg="-use-neg"
use_neg_opt=""
#keep_ann="-keep-ann"
#keep_ann_opt="--keep-ann"
#alg="sinksourceplus"
#goid="GO:0046677"
# best diff Y pestis e0_1
#goid="GO:0016052"
#goid="GO:0071241"
# worst diff for Localplus to SS e-50, Y pestis
#goid="GO:0042221"
# worst diff for Localplus to SS e-25, Y pestis
#goid="GO:0006629"
#goid="GO:0046677"
# best diff 1e-25 M. tuberculosis
#goid="GO:0070566"
#goid="GO:0046777"
#taxon="632"
#taxon="83332"
#taxon="83333"
taxon="208964"
#taxon="99287"
#taxon="243277"
# strange e-coli results
#goid="GO:0006464"
#goid="GO:0006518"
#goid="GO:0009405"
#goid="GO:0006810"
# also look at S. aureus
#taxon="93061"
#goid="GO:0006355"
#goid="GO:0006631"
#goid="GO:0006633"
#declare -a goids=("GO:0072330" "GO:0070566")
#declare -a goids=("GO:1990837")
#declare -a goids=("GO:0006810")  # transport
declare -a goids=("GO:0009432")   # SOS response
#declare -a goids=("GO:0044119" "GO:0044117")
#declare -a version_list=("2018_06-seq-sim-e1e-25" "2018_06-seq-sim-e0_1")
#declare -a version_list=("2018_06-seq-sim-e1e-15" "2018_06-seq-sim-e1e-6")
#declare -a version_list=("2018_06-seq-sim-e0_1") 
declare -a version_list=("2018_06-seq-sim-e0_1-string-700") 
string="-core"
string_opts="--weight-swsn --string$string"
exp_name="${ev_codes}-50-1000${eval_ev_codes}${string}-bp"
#version="2018_06-seq-sim-e1e-25"
#version="2018_06-seq-sim-e0_1"
#opts="--eps 0.0001 --alpha 1.0"
# try using a slightly lower alpha so the propagation stands out a little more
#opts="--eps 0.0001 --alpha 0.99"
opts="--max-iters 10 --eps 0 --alpha 0.95"
#post_opts="--name-postfix -a0_99-2 --num-neighbors 2"
#--name-postfix "-2" \

for goid in ${goids[@]}; do
for version in ${version_list[@]}; do
    cmd="""python src/algorithms/gain_scipy/eval_leave_one_species_out.py \
        --version $version \
        --exp-name $exp_name \
        --pos-neg-file inputs/pos-neg/${ev_codes}/pos-neg-bp-${ann_cutoff}-list.tsv \
        -G $goid \
        -T $taxon --write-prec-rec $use_neg_opt $keep_ann_opt $string_opts \
        -W -1 --forcealg $opts \
        -A localplus -A sinksource -A genemania
    """
        #--pos-neg-file inputs/pos-neg/${ev_codes}/pos-neg-mf-${ann_cutoff}-list.tsv \
        #-A localplus -A sinksource -A genemania
    echo $cmd
    #$cmd

    # now post to graphspace
    #declare -a algs=("localplus" "sinksource" "genemania") 
    #declare -a algs=("localplus") 
    declare -a algs=("sinksource") 
    for alg in ${algs[@]}; do
    cmd="""python2 src/fungcat_post_to_graphspace.py \
        --version $version \
        --net-file inputs/$version/taxon/${taxon}-${exp_name}-swsn.npz
        --exp-name ${exp_name}${use_neg}${keep_ann} \
        --algorithm $alg \
        --goid $goid \
        -U jeffl@vt.edu -P f1fan  \
        --pos-neg-file inputs/pos-neg/${ev_codes}/pos-neg-bp-${ann_cutoff}-list.tsv \
        --pos-neg-file inputs/pos-neg/${ev_codes}/pos-neg-mf-${ann_cutoff}-list.tsv \
        -T $taxon $keep_ann_opt $opts $post_opts \
        --goid-summary-file inputs/pos-neg/${ev_codes}/pos-neg-${ann_cutoff}-summary-stats.tsv
    """
        #--name-postfix "-2" \
        #--num-neighbors 2 \
        #--exp-name ${exp_name}-use-neg \
    echo $cmd
    #$cmd

    done
done
done
