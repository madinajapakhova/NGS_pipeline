#!/bin/bash
analysis_dir=.../analysis_folder/...
bam_files=$analysis_dir/3_Picard_processed/3.3_markDuplicates
filename_extention=.bam

echo 'Started'

# If you need to create a text file with samplename per line
#find $bam_files -type f -name "*.bam" -exec basename {} \; | sed 's/\.[^.]*$//' > $bam_files/bam_files.txt

# Loop over each line in bam_files.txt
for s in $(cat $bam_files/bam_files.txt)
    do
        sleep 1
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.qualimap ./main_child_scripts/child_4_RNASeq_qualimap.sh $analysis_dir $bam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

