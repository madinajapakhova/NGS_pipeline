#!/bin/bash
set -e
analysis_dir=/group/bioinf/Data/Madina_test/pipeline/v1
input_vcf=$analysis_dir/5_VariantCalling/5.3_Genotyping/variants.vcf

echo "Started at `date`"

echo "Started variants filtering"
sbatch ./E.Coli_variant_calling/main_child_scripts/child_5_Filtering.sh $analysis_dir $input_vcf
    
echo "Finished at `date`"
echo "__DONE__"


