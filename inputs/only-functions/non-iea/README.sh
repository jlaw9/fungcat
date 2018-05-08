
# commands I used to get the GO term splits
h="bp"; H="P"; cut1=50; cut2=200; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt
h="mf"; H="F"; cut1=50; cut2=200; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt
h="bp"; H="P"; cut1=200; cut2=500; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt
h="mf"; H="F"; cut1=200; cut2=500; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt
h="bp"; H="P"; cut1=500; cut2=1000; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt
h="mf"; H="F"; cut1=500; cut2=1000; \
    head -n 1 ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv > non-iea-$h-$cut1-$cut2.txt; \
    cat ../../pos-neg/noniea/nonieapos-neg-100-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | grep -P "\t$H\t" | sed "s/GO:0*//" >> non-iea-$h-$cut1-$cut2.txt

