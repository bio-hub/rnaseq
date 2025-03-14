ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

samples=$(cat ./GSE95078.csv | awk -F"," 'NR>1 {print $1}')

for id in ${samples}; do

    echo $id
    prefetch $id
    fasterq-dump --progress --threads $ncores --mem $memory_mega $id
    pigz $id"_1.fastq"
    pigz $id"_2.fastq"
    rm -r $id

done
