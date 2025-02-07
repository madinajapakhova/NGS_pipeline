#!/bin/bash
set -e 
analysis_dir=.../analysis_dir/...
sam_files=$analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/
filename_extention=.sam

echo 'Started'

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/sam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.picard ./main_child_scripts/child_3_Picard.sh $analysis_dir $sam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

