#!/bin/bash
set -e
analysis_dir=/group/bioinf_biomarkers_rna/rna_us_analysis

echo "Started at `date`"

        sbatch -J Annotating_bacteria  ./E.Coli_variant_calling/main_child_scripts/child_Annotating_bacteria.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"
