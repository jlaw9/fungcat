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

The gain-scipy algorithms use the positive and negative examples for each GO term from the [go_term_prediction_examples.py](https://github.com/IGACAT/go_term_prediction_examples) script.

Example call:
```
python src/igacat/go_term_prediction_examples/go_term_prediction_examples.py   \
  --obo-file /data/inputs/goa/2017-09-26-go.obo   \
  --gaf-file inputs/goa/taxon/19-strains-goa.gaf   \
  --out-pref inputs/pos-neg/rem-neg-iea/2018-05-08- \
  --rem-neg-ec IEA --cutoff 50 
```
