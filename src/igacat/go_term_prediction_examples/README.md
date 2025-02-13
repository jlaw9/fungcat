This script takes the annotations in a GAF file, and the GO DAG in an OBO file (available here: http://geneontology.org/page/download-ontology) and assigns every gene as either a positive, negative or unknown for each GO term that passes some criteria. 

For a GO term _t_, we define a gene _g_ as a 
- _positive_ if _g_ is directly annotated to _t_ or to an ancestor of _t_ (more specific term) in the GO DAG
- _negative_ if _g_ is not annotated to _t_ or an ancestor or descendant of _t_ in the GO DAG, but also has at least 1 other annotation
- _unknown_ if _g_ is neither a _positive_ nor a _negative_, meaning it has no annotations, or is annotated to a descendant of _t_ (more general term) in the GO DAG

We say that _g_ is directly annotated to _t_ if the annotation is specified in the GAF file.

The evidence codes used to define the positive, negative and unknown examples can be refined using the 
`--pos-neg-ec`, `--rem-neg-ec` and `--ignore-ec` options described below.

Some examples of how to use these options are in outputs/README.md.


## Usage
```
$ python go_term_prediction_examples.py 
Usage: go_term_prediction_examples.py [options] 
This script takes the annotations in a GAF file, and the GO DAG and assigns
every gene as either a positive (1), negative (-1) or unknown (0) for each GO
term with > cutoff annotations. Writes two tab-separated tables containing the
assignments, one for BP and one for MF, where the rows are genes,  and the
columns are GO term IDs. Also writes a summary statistics table

Options:
  -h, --help            show this help message and exit
  -g GAF_FILE, --gaf-file=GAF_FILE
                        File containing GO annotations in GAF format. Required
  -b OBO_FILE, --obo-file=OBO_FILE
                        GO OBO file which contains the GO DAG. Required
  -c CUTOFF, --cutoff=CUTOFF
                        GO terms having > cutoff positive instances (proteins)
                        are kept. Default=1000
  -o OUT_PREF, --out-pref=OUT_PREF
                        Prefix used to write a table of positives, negatives,
                        and unknowns for each GO category.Writes an output
                        file for BP and MF: <out-pref>pos-neg-<cutoff>-P.tsv
                        and <out-pref>pos-neg-<cutoff>-F.tsv
  --pos-neg-ec=POS_NEG_EC
                        Comma-separated list of evidence codes used to assign
                        positive and negative examples. If none are specified,
                        all codes not in the two other categories (--rem-neg-
                        ec and --ignore-ec) will be used by default.
  --rem-neg-ec=REM_NEG_EC
                        Comma-separated list of evidence codes used to remove
                        negative examples. Specifically, If a protein would be
                        labelled as a negative example for a given term but is
                        annotated with a 'rem_neg' evidence code for the term,
                        it is instead labelled as unknown. If none are
                        specified, but --pos-neg-ec codes are given, all codes
                        not in the other two categories will be put in this
                        category by default.
  --ignore-ec=IGNORE_EC
                        Comma-separated list of evidence codes where
                        annotations with the specified codes will be ignored
                        when parsing the GAF file. For example, specifying
                        'IEA' will skip all annotations with an evidence code
                        'IEA'. If both --pos-neg-ec and --rem-neg-ec codes are
                        given, everything else will be ignored by default.

--gaf-file (-g), --obo-file (-b), and --out-pref (-o) are required
```

## Installation
This script requires Python 3 due to the use of obonet to build the GO DAG.

Required packages: obonet, networkx, pandas, tqdm

To install the required packages:
```
pip install -r requirements.txt
```

Optional: use a virtual environment
```
virtualenv -p /usr/bin/python3 py3env
source py3env/bin/activate
pip install -r requirements.txt
```
