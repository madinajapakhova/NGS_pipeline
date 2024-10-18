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
# In this script: 1.3 QUALITY CONTROL of the CURATED DATA: generate FastQC and MultiQC reports
#         - run fastqc on each .fastq file
#         - run multiqc on all fastqc reports 
#         - if paired-end, compare quality of forward and reverse reads  
# ******************************************************************************


# ******************************************************************************
#                           !USER PROVIDED HEADING!
# ******************************************************************************
analysis_dir=.../analysis_directory_example
paired_end="Yes"  
fastQ_files_folder=.../analysis_directory_example/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters
fastQ_files_R1=$fastQ_files_folder/R1
fastQ_files_R2=$fastQ_files_folder/R2 
sym_links="No"
extension='*fastq.gz'
sw=.../FastQC_sw_example/FastQC
path_to_java=.../path_to_java_example
threads=4 
memory_per_file=512 

# analysis_dir --- absolute path to the folder where each step of pipeline execution is stored. Each step is stored in a separate subfolder.
# paired_end --- "Yes" if you want to additionally compare the quality of forward and reverse reads.
# Input: path to .fastq files 
# !filenames should not contain spaces!
# fastQ_files_folder --- where the curated .fastq data is kept. 
# Paired-end case: forward and reverse reads in separate subfolders - 'R1' and 'R2'
# fastQ_files_R1 --- for paired-end reads, subdirectory for the forward reads. Has path $fastQ_files_folder/R1 
# fastQ_files_R2 --- for paired-end reads, subdirectory for the reverse reads. Has path $fastQ_files_folder/R2
# sym_links --- if you create links for the .fastq.gz files and provide these symlinks as input, then set this argument to "Yes"
# extension --- extension of your fastq files
# threads --- # number of files which can be processed simultaneously.  Each thread will be allocated 250MB of memory, so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine
# memory_per_file --- Base amount of memory, in Megabytes, to process each file.
# path_to_java --- Path where you store your java distribution
# ******************************************************************************

echo "Started at `date`"
# ******************************************************************************
# Create a list of files to be passed to FastQC
# !Filenames should not contain spaces!
# ******************************************************************************
if test $(find "$fastQ_files_folder" -name "* *" | wc -c) -gt 0; then
    echo "Error: filenames contain spaces"
    exit 
else
    echo "Check passed. Filenames do not contain spaces"    
fi
    
# create a vector with all samples
if [ "$sym_links" = "Yes" ]; then
    echo "Input files were provided as symlinks"
    files=$(find -L "$fastQ_files_folder" -type f -name "$extension" -print0 | tr '\0' ' ')
else
    files=$(find "$fastQ_files_folder" -type f -name "$extension" -print0 | tr '\0' ' ')
fi

# ******************************************************************************
# Running FastQC
# Installation: https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt
# Syntax: https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
# Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# ******************************************************************************

mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/ -p
mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ -p 

echo "Running FastQC"

$sw/fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ \
    --extract \
    --memory $memory_per_file \
    --threads $threads \
    --java $path_to_java \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_log.txt 2>&1
        
# ******************************************************************************
# Runing multiqc on fastqc reports
# Multiqc aggregates multiple fastqc reports into a single .html report
# Inputs for multiqc are fastqc reports
# https://multiqc.info/docs/
# ******************************************************************************
mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ -p
# -m fastqc: because we're only reporting on FastQC, we specify this module to speed up multiqc
# --interactive: makes the .html file interactive
# -f: overwrites an existing report
echo "MultiQC: aggregating FastQC reports for all samples"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ -o $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_log.txt 2>&1
echo "MultiQC report is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ 

#*******************************************************************************
# In case the reads are paired-end, you might want to run quality checks only on forward reads
# or only on reverse reads
#*******************************************************************************
if [ "$paired_end" = "Yes" ]; then
    echo "Running quality control separately on forward and reverse reads"
    
# forward reads
    echo "First checking only the forward reads"
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ -p
    
    if [ "$sym_links" = "Yes" ]; then
        echo "Input files were provided as symlinks"
        files=$(find -L "$fastQ_files_R1" -type f -name $extension -print0 | tr '\0' ' ')
    else
        files=$(find "$fastQ_files_R1" -type f -name $extension -print0 | tr '\0' ' ')
    fi
    
$sw/fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ \
    --extract \
    --memory $memory_per_file \
    --threads $threads \
    --java $path_to_java \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_fwd_reads_log.txt 2>&1
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/ -p

echo "MultiQC: aggregating FastQC reports for the forward reads"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ -o  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_fwd_log.txt 2>&1
echo "MultiQC report for the forward reads (R1) is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/

# reverse reads
    echo "Checking only the reverse reads"
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads/ -p
    if [ "$sym_links" = "Yes" ]; then
        echo "Input files were provided as symlinks"
        files=$(find -L "$fastQ_files_R2" -type f -name $extension -print0 | tr '\0' ' ')
    else
        files=$(find "$fastQ_files_R2" -type f -name $extension -print0 | tr '\0' ' ')
    fi
      
$sw/fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads/ \
    --extract \
    --memory $memory_per_file \
    --threads $threads \
    --java $path_to_java \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_rev_reads_log.txt 2>&1
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/ -p

echo "MultiQC: aggregating FastQC reports for the reverse reads"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads/ -o $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_rev_log.txt 2>&1
echo "MultiQC report for the reverse reads (R2) is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/

else
    echo "Additional quality checks are not needed"
fi

echo "Finished at `date`"
echo "__DONE__"

    
    
    
    
    
    
    
    
