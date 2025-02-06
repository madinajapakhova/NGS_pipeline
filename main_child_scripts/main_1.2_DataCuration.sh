#!/bin/bash
set -e

# ******************************************************************************
#   !USER-PROVIDED INPUT. LEAVE NONE EMPTY!
# ******************************************************************************

analysis_dir=.../analysis_folder/...
folder_nextstep=.../raw_data_folder/...
paired_end="Yes"  
filename_list_exclude=.../raw_data_folder/exclude_list.txt 
softwares=.../software_folder
UMI_extraction="No"  
UMI_pattern="None" 
UMI_pattern2="None" 
umi_method=string 
cut_adapters="Yes" 
trimming_software="cutadapt" 
adapter_fwd="GTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTA" 
adapter_rev="GTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTA" 
illumina_adapters_file_PE="None" 
illumina_adapters_file_SE="None"
low_quality="No" 
remove_poly="No" 
threads=4 
readname=_R
filename_extention=_001.fastq.gz 

              
# analysis_dir --- absolute path to the folder where each step of pipeline execution is stored. Each step is stored in a separate subfolder.
# folder_nextstep --- absolute path for storing most recent FASTQ data. Updated after each manipulation on the data
# folder_nextstep is where the most recent .fastq.gz files will be stored. Includes subfolders "R1" and "R2" if reads are paired-end.
# folder_nextstep will be updated as we move from one process to another; if process is skipped folder_nextstep is not updated
# starting folder_nextstep is where our .fastq.gz files ended up in 1.1_QualityCheckRaw: $analysis_dir/0_fastQ_data
# paired_end --- "Yes" for paired end reads; "No" otherwise 
# filename_list_exclude --- a .txt file with sample names that need to be excluded 
# filename_extension --- extension of the file to work on. Data curation works on .fastq.gz
# samplename --- name of the sample without extension 
# readname --- read naming convention, e.g. R
# filename_list_exclude --- if you want to exclude some samples from the pipeline provide a list as a text file. One sample per line
# softwares --- folder with tools
# UMI_extraction --- "Yes" if UMIs apply; "No" otherwise 
# UMI_pattern --- UMIs pattern for forward reads if applies, e.g. "NNN". Set to "None" if does not apply   
# UMI_pattern2 --- UMIs pattern for reverse reads if applies, e.g. "NNN". Set to "None" if does not apply
# umi_method --- UMIs extraction method "string" or "regex". Set to "None" if does not apply
# cut_adapters --- "Yes" if removing adapters applies; "No" otherwise 
# trimming_software --- "cutadapt" or "trimmomatic" depending on which tool you prefer for your data. Trimmomatic is generally better for paired-end Illumina reads  
# adapter_fwd --- adapter sequence for forward reads  Set to "None" if does not apply
# adapter_rev --- adapter sequence for reverse reads. Set to "None" if does not apply
# illumina_adapters_file_PE --- for trimmomatic a file with adapters needs to be passed. PE = paired end. Set to "None" if does not apply
# illumina_adapters_file_SE --- for trimmomatic a file with adapters needs to be passed. PE = single end. Set to "None" if does not apply
# low_quality --- "Yes" if trimming low quality reads applies; "No" otherwise
# remove_poly --- "Yes" if removing poly-sequences applies; "No" otherwise
# poly_seq --- poly-sequences to be removed, e.g. for poly_A set to "A{50}", for wildcard set to "N{50}". Set to "None" if does not apply
# threads --- some softwares, e.g cutadapt and trimmomatic support multithreading, i.e. when one samples is processed simultaneously on multiple threads. Specify number of threads. 

echo "Started at `date`"
# ******************************************************************************
# Sanity check
# ******************************************************************************
echo 'Checking that all params have been provided'
params=("$analysis_dir" "$folder_nextstep" "$paired_end" "$filename_list_exclude" "$softwares" "$UMI_extraction" "$UMI_pattern" "$UMI_pattern2" "$umi_method" "$cut_adapters" "$trimming_software" "$adapter_fwd" "$adapter_rev" "$illumina_adapters_file_PE" "$illumina_adapters_file_SE" "$low_quality" "$remove_poly" "$threads" "$readname" "$filename_extention")

for param in "${params[@]}"; do
  if [ -z "$param" ]; then
    echo "Error: One or more parameters are empty."
    exit 
  else echo "Input has been provided. All good."
  fi
done

# ******************************************************************************
# Provide a text file with samplenames (without extension). One sample per line
# ******************************************************************************
# Each sample has to follow the following naming convention: samplename + readname + extension
# For single-end reads, provide the full readname, e.g. R1
# For paired-end reads, provide the readname without 1 and 2. The pipeline automatically appends 1 to forward reads, and 2 to reverse. 

# Example 1. You have a sample with paired-end reads "DNA1_R1.fastq.gz" and "DNA1_R2.fastq.gz"
# Then, samplename is DNA1; readname is _R; extension is .fastq.gz 
# The pipeline appends automatically 1 and 2 to the filename 
# So the result will be two constructed filenames: DNA1_R1.fastq.gz and DNA1_R2.fastq.gz

# Example 2. Sample with paired-end reads "DNA1_R1_001.fastq.gz" and "DNA1_R2_001.fastq.gz" 
# Samplename is DNA1, readname is _R, and extension would be _001.fastq.gz
# Constructed filenames are then DNA1_R1_001.fastq.gz and DNA1_R2_001.fastq.gz 

# Example 3. Sample with single-end reads "DNA1_R1.fastq.gz" 
# samplename - DNA1, readname - _R1

# ******************************************************************************
# Activate the child script
# ******************************************************************************

#chmod u+x ./main_child_scripts/child_1.2_DataCuration.sh

# this can help you create a text file with a list if samplenames. 
# Beware! This is highly specific to the filename conventions used in your dataset. So take this line more as a hint!
# find . -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/_R.*//' | sort | uniq > $analysis_dir/samplenames.txt

# ******************************************************************************
# Run the pipeline for each sample in the list
# ******************************************************************************
for s in $(cat $analysis_dir/samplenames.txt)
    do
        sleep 1
        samplename="$s"       
        echo "Processing sample: $samplename"
        sbatch -J 1.2_DataCuration.$samplename ./main_child_scripts/child_1.2_DataCuration.sh $analysis_dir $folder_nextstep $paired_end $filename_list_exclude $softwares $UMI_extraction $UMI_pattern $UMI_pattern2 $umi_method $cut_adapters $trimming_software $adapter_fwd $adapter_rev $illumina_adapters_file_PE $illumina_adapters_file_SE $low_quality $remove_poly $threads $readname $filename_extention  $samplename 
    done

echo "Finished at `date`"
echo "__DONE__"
