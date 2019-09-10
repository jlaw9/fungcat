# for each species,
# store the SWSN weighted network, then compute the distsances
source virtual-envs/py3env/bin/activate


#taxon="83332"; 
#taxon="83333"; 
taxon="208964"; 
version="2018_06-seq-sim-e0_1-string-700";
# expc
# first create the SWSN network
# don't need to actually run sinksource, but I don't have a way to only generate the weighted network currently
#python src/algorithms/gain_scipy/eval_leave_one_species_out.py \
#	--version $version --exp-name weight-test \
#	--pos-neg-file inputs/pos-neg/expc/pos-neg-bp-10-list.tsv  \
#	--only-functions inputs/only-functions/expc/expc-50-1000.txt  -T $taxon \
#	--weight-swsn --string-core -W 0  \
#	--max-iters 10 --eps 0 --alpha 0.95 -A sinksource  \
#	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-50-1000-core-bp-swsn.npz \
#	--num-test-cutoff 10 --forcealg
# then compute the weighted shortest paths
python src/analyze_results/analyze_leave_one_speces_out.py \
	--version $version --exp-name expc-50-1000-core-bp \
	--pos-neg-file inputs/pos-neg/expc/pos-neg-bp-10-list.tsv  \
	--only-functions inputs/only-functions/expc/expc-50-1000.txt  -T $taxon \
	--string-core -W 0  \
	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-50-1000-core-bp-swsn.npz \
	--num-test-cutoff 10

# COMP
#python src/algorithms/gain_scipy/eval_leave_one_species_out.py \
#	--version $version --exp-name weight-test \
#	--pos-neg-file inputs/pos-neg/expc-comp/pos-neg-bp-10-list.tsv \
#	--pos-neg-file-eval inputs/pos-neg/comp/pos-neg-bp-5-list.tsv \
#	--only-functions inputs/only-functions/expc-comp/expc-comp-50-1000.txt  -T $taxon \
#	--weight-swsn --string-core -W 0  \
#	--max-iters 10 --eps 0 --alpha 0.95 -A sinksource  \
#	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-comp-50-1000comp-core-bp-swsn.npz \
#	--num-test-cutoff 10 --forcealg
python src/analyze_results/analyze_leave_one_speces_out.py \
	--version $version \
	--exp-name expc-comp-50-1000comp-core-bp \
	--pos-neg-file inputs/pos-neg/expc-comp/pos-neg-bp-10-list.tsv \
	--pos-neg-file-eval inputs/pos-neg/comp/pos-neg-bp-5-list.tsv \
	--only-functions inputs/only-functions/expc-comp/expc-comp-50-1000.txt  -T $taxon \
	--string-core -W 0  \
	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-comp-50-1000comp-core-bp-swsn.npz \
	--num-test-cutoff 10

# ELEC
#python src/algorithms/gain_scipy/eval_leave_one_species_out.py \
#	--version $version \
#	--exp-name weight-test \
#	--pos-neg-file inputs/pos-neg/expc-comp/pos-neg-bp-10-list.tsv  \
#	--pos-neg-file-eval inputs/pos-neg/iea/pos-neg-bp-10-list.tsv \
#	--only-functions inputs/only-functions/expc-comp/expc-comp-50-1000.txt  -T $taxon \
#	--weight-swsn --string-core -W 0  \
#	--max-iters 10 --eps 0 --alpha 0.95 -A sinksource  \
#	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-comp-50-1000iea-core-bp-swsn.npz \
#	--num-test-cutoff 10 --forcealg
#python src/analyze_results/analyze_leave_one_speces_out.py \
#	--version $version \
#	--exp-name expc-comp-50-1000iea-core-bp \
#	--pos-neg-file inputs/pos-neg/expc-comp/pos-neg-bp-10-list.tsv \
#	--pos-neg-file-eval inputs/pos-neg/iea/pos-neg-bp-10-list.tsv \
#	--only-functions inputs/only-functions/expc-comp/expc-comp-50-1000.txt  -T $taxon \
#	--string-core -W 0  \
#	--net-file inputs/2018_06-seq-sim-e0_1-string-700/taxon/$taxon-expc-comp-50-1000iea-core-bp-swsn.npz \
#	--num-test-cutoff 10

