for folder in neanderthals_four
do
dir=/data/coxvgi/zhoud2/projects/cross_tissue/v6p/out
mkdir ${dir}/result/${folder}
i_dir=${dir}/${folder}
rm ${dir}/tmp/*
comb_list=$(ls ${i_dir})
for comb in $comb_list
do
echo $comb
cat ${i_dir}/${comb}/* > ${dir}/tmp/${comb}.txt
head -n 1 ${i_dir}/${comb}/1_* > ${dir}/tmp/head.txt
cat ${dir}/tmp/${comb}.txt | awk -F"\t" '$1!="geneid"' > ${dir}/tmp/${comb}_nohead.txt
cat ${dir}/tmp/head.txt ${dir}/tmp/${comb}_nohead.txt > ${dir}/result/${folder}/${comb}
done
done
