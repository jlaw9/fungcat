# fungcat

## Setup on the bioinformatics filesystem
The GO Annotations (GOA) files are located at `/data/inputs/goa`
The STRING networks are located at `/data/inputs/string`
To use those files, add a symbolic link from your inputs folder to those locations:
```
cd inputs
ln -s /data/inputs/goa goa
ln -s /data/inputs/string string
```

The Biorithm/GAIN package can be downloaded from https://bioinformatics.cs.vt.edu/~murali/software/biorithm/index.html

Calling the `master-script.py` will setup most of the inputs for you. If you specify an `--exp-name`, it will try to run my copy of GAIN with the default settings.

Example call:
```
python master-script.py --version 2017_10-seq-sim
```

The gain-scipy algorithms use the positive and negative examples for each GO term from the [go_term_prediction_examples.py](https://github.com/IGACAT/go_term_prediction_examples) script.
TODO add this to `master-scripty.py`.

Example call:
```
python src/igacat/go_term_prediction_examples/go_term_prediction_examples.py   \
  --obo-file /data/inputs/goa/2017-09-26-go.obo   \
  --gaf-file inputs/goa/taxon/19-strains-goa.gaf   \
  --out-pref inputs/pos-neg/rem-neg-iea/2018-05-08- \
  --rem-neg-ec IEA --cutoff 50 
```

From there, you can run the gain-scipy scripts using a call like this:
```
python src/algorithms/gain-scipy/run_algs.py \
	--version 2017_10-seq-sim  \
	--exp-name rem-neg-iea-40-50 \
	--pos-neg-file inputs/pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-bp-10-list.tsv \
	--pos-neg-file inputs/pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-mf-10-list.tsv \
	--only-functions inputs/only-functions/rem-neg-iea/rem-neg-iea-40-50.txt -A sinksourceplus-ripple -A sinksourceplus-squeeze  -k 200 -a 0.8 \
	--epsUB 0.01 \
	--num-pred-to-write 200 \
	--forcealg \
	--verbose
```
