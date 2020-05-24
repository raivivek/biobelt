#! /bin/bash

help="Usage: gregorRun.sh <output_dir> <bed_files>.."

[[ -z $1 ]] && { echo $help; exit 1; }
[[ -z $2 ]] && { echo $help; exit 1; }

out_dir="$1"

## create index bed
index_bed=$out_dir/indexbed.txt

# GWAS-SNP files
echo "drmr:job processors=6 job_name=gregor"
for i in ~vivekrai/data/gwas_arushi/*.txt; do
  name=$(basename $i | sed 's/.txt//g')

  echo "mkdir -p $out_dir/traits/$name && \
    gregorMakeConfig -d $out_dir/traits/$name/output \
    --annotfiles "$2" \
    --indexbed $index_bed \
    --outfile $out_dir/traits/$name/$name.conf \
    --snpfile $i && ionice -c2 -n7 GREGOR.pl --conf $out_dir/traits/$name/$name.conf"
done
