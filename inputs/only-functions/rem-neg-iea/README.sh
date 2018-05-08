
# use this command to get results for all goterms
cut1=30; cut2=40;
head -n 1 ../../pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-10-summary-stats.tsv > rem-neg-iea-$cut1-$cut2.txt;
sed -i "s/^/#/" rem-neg-iea-$cut1-$cut2.txt;
cat ../../pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-10-summary-stats.tsv | awk -v cut1=$cut1 -v cut2=$cut2 -F'\t' '{if($4 >= cut1 && $4 < cut2){ print $0}}' | sed "s/GO:0*//" >> rem-neg-iea-$cut1-$cut2.txt
