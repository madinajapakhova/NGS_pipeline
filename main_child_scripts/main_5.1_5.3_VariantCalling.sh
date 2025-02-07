#!/bin/bash
set -e
analysis_dir=.../analysis_folder/... 
bam_files=$analysis_dir/3_Picard_processed/3.3_markDuplicates                                        
filename_extention=.bam

echo 'Started'

# Loop over each line in samplenames.txt

for s in $(cat $bam_files/bam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J gatk_gvcf.$samplename ./main_child_scripts/child_5.1_5.3_VariantCalling.sh $analysis_dir $bam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"
