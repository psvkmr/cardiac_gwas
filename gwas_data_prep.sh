for i in `seq 1 22`
do
  sort -k2 -g GWAS*Chr${i}_2*.linear > chr${i}_sumstats_sorted.txt
done

ls -v chr*sorted.txt | xargs cat > total_sumstats_sorted.txt
