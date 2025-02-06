#!/bin/bash
set -e

analysis_dir=.../analysis_directory_example
processed_files=$analysis_dir/curated_data_directory                                      
paired_end="Yes" 
filename_extention=_001.fastq.gz 
readname=_R
genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir
reference_genome38=.../resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta #FASTA file with reference genome (mandatory)
reference_gtf=.../resources/gencode.v38.annotation.gtf #GTF file with annptated transcripts (highly recommended)
threads=4
overhang=100

echo 'Started'

# ******************************************************************************
# Sanity check
# ******************************************************************************
echo 'Checking that all params have been provided'
params=("$analysis_dir" "$processed_files" "$paired_end" "$filename_extention" "$readname" "$genomeDir" "$reference_genome38" "$reference_gtf" "$threads" "$overhang")

for param in "${params[@]}"; do
  if [ -z "$param" ]; then
    echo "Error: One or more parameters are empty."
    exit 
  else echo "Input has been provided. All good."
  fi
done

# ******************************************************************************
# Activate the child script
# chmod u+x ./child/child_1.2_DataCuration.sh
# ******************************************************************************

# this can help you create file with samplenames. 
# Beware! This is highly specific to the filename conventions used in your dataset. So take this line more as a hint!
# find . -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/_R.*//' | sort | uniq > $analysis_dir/samplenames.txt
   
echo "2.1. Generating reference genome"
for s in $(cat $analysis_dir/samplenames.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.2_mapping ./main_child_scripts/child_2_Mapping_RNA.sh $analysis_dir $processed_files $paired_end $filename_extention $samplename $readname $genomeDir $reference_genome38 $reference_gtf $threads $overhang
    done

echo "Finished at `date`"
echo "__DONE__"
