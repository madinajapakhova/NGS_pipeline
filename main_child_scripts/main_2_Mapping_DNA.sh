#!/bin/bash
#set -x
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/bacteria
fastQ_files_folder=/group/bioinf/Users/Madina/pipeline_tested/bacteria/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed                          
paired_end="Yes" 
filename_extention=.fastq.gz
readname=R

echo 'Started'

#!/bin/bash

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/samplenames.txt)
    do
	sleep 1
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J 2_Align.$samplename ./main_child_scripts/child_2_Mapping_DNA.sh $analysis_dir $fastQ_files_folder $paired_end $filename_extention $samplename $readname
    done

echo "Finished at `date`"
echo "__DONE__"
