This directory contains various sets of GO terms used to inform GAIN which GO terms
to generate predictions/evaluate. The leading "GO:00" of GO IDs is removed as
GAIN doesn't recognize GO terms with them (see below).

Dec 21 2017:
These files contain the GO terms with > 1000 proteins annotated to them. Copied from inputs/pos-neg/19-strainspos-neg-1000-summary-stats.tsv and inputs/pos-neg/non-iea-pos-neg-1000-summary-stats.tsv
1000-annotations.txt
non-iea-1000-annotations.txt

Then I removed the GO: with leading 0s to the IDs because GAIN won't recognize
the GO terms with them
sed -i "s/GO:0*//g" 1000-annotations.txt
sed -i "s/GO:0*//g" non-iea-1000-annotations.txt

