
## Setup sparse matrices
Example call to setup sparse matrices for the core STRING networks and the GO annotations for 50-1000 proteins annotated with non-IEA.

```sh
python src/algorithms/setup_sparse_networks.py  \
  --version 2017_10-string \
  --pos-neg-file inputs/pos-neg/rem-neg-iea/pos-neg-bp-50-list.tsv  \
  --only-functions inputs/only-functions/rem-neg-iea/rem-neg-iea-50-1000.txt \
  --out-pref-net inputs/2017_10-string/sparse_nets/ \
  --out-pref-ann inputs/2017_10-string/rem-neg-iea/bp-50-1000- \
  --core-string
```

## Weight networks using SWSN
To combine the networks using the SWSN (Simultaneous Weighting with Specific Negatives) method,
use the `--weight-SWSN` option, and then pass the resulting netork to the 
`run_algs.py` script for cross-validation.

```sh
python src/algorithms/setup_sparse_networks.py  \
  --version 2017_10-seq-sim-string-swsn \
  --pos-neg-file inputs/pos-neg/rem-neg-iea/pos-neg-bp-50-list.tsv  \
  --only-functions inputs/only-functions/rem-neg-iea/rem-neg-iea-50-1000.txt \
  --out-pref-net inputs/2017_10-seq-sim-string-swsn/sparse_nets/ \
  --out-pref-ann inputs/2017_10-seq-sim-string-swsn/rem-neg-iea/bp-50-1000- \
  --only-core \
  --weight-SWSN inputs/2017_10-seq-sim-string-swsn/rem-neg-iea/bp-50-1000-
```

## Run cross-validation using the combined network
```sh
python -u src/algorithms/gain_scipy/run_algs.py \
  --version 2017_10-seq-sim-string-swsn \
  --net-file inputs/2017_10-seq-sim-string-swsn/rem-neg-iea/50-1000-7-nets-combined-SWSN.npz  \
  --exp-name rem-neg-iea-bp-50-1000-core \
  --pos-neg-file inputs/pos-neg/rem-neg-iea/pos-neg-bp-50-list.tsv \
  --pos-neg-file inputs/pos-neg/rem-neg-iea/pos-neg-mf-50-list.tsv \
  --only-functions inputs/only-functions/rem-neg-iea/rem-neg-iea-50-1000.txt  \
  -A sinksource -A localplus \
  -a 0.8 --eps 0 --max-iters 20 \
  --only-cv -C 5
```
