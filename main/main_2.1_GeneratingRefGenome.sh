#!/bin/bash
set -e

# ******************************************************************************
#   !USER-PROVIDED INPUT. LEAVE NONE EMPTY!
# ******************************************************************************
analysis_dir=.../analysis_directory_example
threads=4
reference_genome38=.../resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 
reference_gtf=.../gencode.v38.annotation.gtf 
temp_folder=.../star 

# analysis_dir --- absolute path to the folder where each step of pipeline execution is stored. Each step is stored in a separate subfolder.
# reference_genome38 --- FASTA file with reference genome (mandatory)
# reference_gtf --- GTF file with annptated transcripts (highly recommended)
#   !MAKE SURE YOUR FASTA AND GTF FILES ARE COMPATIBLE, THAT IS, THEY USE SAME CHROMOSOME NAMING CONVENTIONS (CHR1, CHROM1, 1, ETC)!        
#   For reference genome, this scripts uses hg38 from the GATK bundle. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
#   For GTF file, it used the gencode gtf file.  See https://www.gencodegenes.org/human/release_38.html 

echo "Started at `date`"

# ******************************************************************************
# Sanity check
# ******************************************************************************
echo 'Checking that all params have been provided'
params=("$analysis_dir" "$threads" "$reference_genome38" "$reference_gtf" "$temp_folder")

for param in "${params[@]}"; do
  if [ -z "$param" ]; then
    echo "Error: One or more parameters are empty."
    exit 
  else echo "Input has been provided. All good."
  fi
done

      
echo "2.1. Generating reference genome"
sbatch -J 2.1_GeneratingRefGenome $analysis_dir/child/child_2.1_GeneratingRefGenome.sh $analysis_dir $threads $reference_genome38 $reference_gtf $temp_folder
   
echo "Finished at `date`"
echo "__DONE__"

