#!/usr/bin/env bash
set -e 

# ******************************************************************************
#                   PIPELINE OVERVIEW
# ******************************************************************************
# I DATA CLEANING (raw FASTQ -> curated FASTQ)
# 1.1 Quality check of the raw data
# 1.2 Data Curation. Depending on the needs of the data:
#        - excluding samples 
#        - UMI extraction 
#        - cutting adapter sequences 
#        - remove low quality reads 
#        - removing poly-sequences 
# 1.3 Quality check of the curated data

# II Alignment -> (curated FASTQ -> SAM)
# 2.1 Generating genome indexes files. (Needs to be done only once!) 
# 2.2 - 2.4 Alignment 

# III Picard Processing (SAM -> BAM)
# 3.1 Sorting
# 3.2 Adding read group information
# 3.3 Marking duplicates and indexing

# IV Quality check of aligned bam files
# 4.1 Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: 1) sample was contaminated (-> check .fastq data quality) 2) alignment did not work (-> check mapping options)

# 4.2 Qualimap
#   BAM quality check
#   RNA quality check if working with RNA-Seq

# V Variant Calling (BAM -> VCF)
# 5.1 SplitNCigarReads
# 5.2.a BaseRecalibrator
# 5.2.b ApplyBQSR
# 5.3 HaplotypeCaller
# 5.4 GenomicsDBImport
# 5.5 GenotypeGVCFs
# 5.6 Merging per/chromosome VCFs into one
# 5.7 VQSR  
# 5.8 Hard filtering
# ******************************************************************************

# ******************************************************************************
# In this script: # 2.1 Generating genome indices files
#                       RNA-STAR version 2.7.9a
#                       Reference genome version: GRCh38.p13, see https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/ (RefSeq)
# ******************************************************************************

# ******************************************************************************
#                           !USER DEFINED HEADING!
# ******************************************************************************
# External shell arguments:
# ******************************************************************************
analysis_dir=$1 
threads=$2
reference_genome38=$3
reference_gtf=$4 
temp_folder=$5
# ******************************************************************************

if [ -d "$temp_folder" ]; then
  rm -r $temp_folder
fi

#*******************************************************************************
# 2.1 Generating genome indexes files
#*******************************************************************************

# Important: generating indices files and mapping have to be done with same version of STAR

# Directory where the indices files will be stored. For human genome we need to have at least 100GB space on disk.
mkdir $analysis_dir/2_mapping_STAR/ -p
mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir/ -p
genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir


echo "Starting step 2.1 Generating genome indexes files"

#STAR manual https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf

#   --genomeChrBinNbits: if working with a large (>5000 references), this helps to reduce RAM consumption
#                        recommended min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])

#$softwares/STAR --runMode genomeGenerate \
 STAR --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles  $reference_genome38 \
    --sjdbGTFfile $reference_gtf \
    --sjdbOverhang 75 \
    --runThreadN $threads \
    --outTmpDir $temp_folder


echo "Finished step 2.1 Generating genome indexes files"

