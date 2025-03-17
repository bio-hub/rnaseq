mkdir samples
mkdir stratookit_cache

ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

cd samples

samples=$(cat ./GSE206364.csv | awk -F"," 'NR>1 {print $1}')

for id in ${samples}; do

    echo $id
    prefetch $id
    fasterq-dump --progress --threads $ncores --mem $memory_mega $id
    pigz $id".fastq"

done
