The output in the two folders here `rem-neg-iea` and `exp` were created using the two commands below.

1) Example call to get positive and negative examples where proteins 
that would be labelled as "negatives" but are annotated with an IEA evidence code
to a given term (or descendant of the term) are instead labelled as "unknowns"

``` sh
mkdir outputs/rem-neg-iea
python go_term_prediction_examples.py \
  --gaf-file inputs/19-strains-goa.gaf \
  --obo-file inputs/2017-09-26-go.obo \
  --cutoff 1000 \
  --rem-neg-ec IEA \
  --out-pref outputs/rem-neg-iea/rem-neg-iea- \
  > outputs/rem-neg-iea/2018-01-15-rem-neg-iea-1000.log
```

2) Example call to assign positive and negative examples using experimental evidence codes
where proteins that would be labelled as "negatives" but are annotated with some other evidence code
to a given term (or descendant of the term) are instead labelled as "unknowns"

``` sh
mkdir outputs/exp
python go_term_prediction_examples.py \
  --gaf-file inputs/19-strains-goa.gaf \
  --obo-file inputs/2017-09-26-go.obo \
  --cutoff 250 \
  --pos-neg-ec EXP,IDA,IPI,IMP,IGI,IEP \
  --out-pref outputs/exp/exp- \
  > outputs/exp/2018-01-15-exp-250.log
```
