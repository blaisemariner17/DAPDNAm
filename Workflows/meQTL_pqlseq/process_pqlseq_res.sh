
for i in {1..38}
do
	awk 'FNR>1' pqlseq_results/pqlseq_res${i}_* > processed_files/all_pqlseq_res${i}.tsv
done

cat processed_files/all_pqlseq_res*.tsv > all_pqlseq_res.tsv
rm processed_files/*
