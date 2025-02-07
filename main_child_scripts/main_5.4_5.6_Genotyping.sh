#!/bin/bash
set -e
analysis_dir=.../analysis_folder/...

echo "Started at `date`"
for c in $(cat $analysis_dir/chromosomes.txt) 
    do
        chrom="$c"
        echo "Processing chromosome: $chrom"
        sbatch -J $chrom.Genotyping ./main_child_scripts/child_5.4_5.6_Genotyping.sh $analysis_dir $chrom
    done


echo "Finished at `date`"
echo "__DONE__"


