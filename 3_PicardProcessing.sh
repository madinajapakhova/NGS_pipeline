#!/usr/bin/env bash
set -e 

# ******************************************************************************
#                   OVERVIEW
# ******************************************************************************
# I DATA CLEANING
# 1.1 Quality check of the raw data. FastQC and MultiQC
# 1.2 Data Curation
#        - excluding samples (optional)
#        - UMI extraction (optional; usually NOT needed for bacteria)
#        - cutting adapter sequences (optional)
#        - remove low quality reads (optional)
#        - removing poly-sequences (optional)
# 1.3 Quality check of the curated data

# II Alignment -> .sam file
# 2.1 Generating genome indexes files. (Needs to be done only once!) 
# 2.2 - 2.4 Alignment 

# III Picard Processing -> .bam file
# 3.1 Sorting
# 3.2 Adding read group information
# 3.3 Marking duplicates and indexing

# IV Quality check of aligned bam files
# 4.1 Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: 1) sample was contaminated (-> check .fastq data quality) 2) alignment did not work (-> check mapping options)

# 4.2 Qualimap

# V Variant Calling
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
# In this script: STEP 3. Picard Processing
#       - Adding read group information
#       - Sorting
#       - Marking duplicates and indexing
# ******************************************************************************

# ******************************************************************************
#                           !USER PROVIDED HEADING!
# ******************************************************************************

# External shell arguments:
# ******************************************************************************
analysis_dir=$1 
sam_files=$2
filename_extention=$3 #.sam
samplename=$4 
# ******************************************************************************

# Constructing the name of the file 
filename="${samplename}Aligned.out${filename_extention}"

# Internal arguments:
# ******************************************************************************
sw=.../software_folder/...
memory=-Xmx30g
# ******************************************************************************
# ******************************************************************************

# Picard tools require Java 1.17 or newer
# singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java
# -Xmx60g indicating the amount of memory to allocate to the Java virtual machine which will run the program. 60 gigabytes

# analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders at each each step.
# sw --- folder with java based softwares, e.g. fastqc, picard, cutadapt
# sam_files -- path to where the aligned sam files were saved in the Mapping step
# ******************************************************************************

# ******************************************************************************
# LATEST VERSION PICARD UPDATES: https://github.com/broadinstitute/picard/wiki
#       - PICARD installation 
#           - Picard 3.0 requires Java 1.17 version
#           - Picard is now built using gradle

#       -PICARD syntax
#       -NOTE: Picard's command line syntax has changed
#       -wherever you used to do e.g. I=input.bam, you now do -I input.bam.
#       -see https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
# ******************************************************************************

mkdir $analysis_dir/3_Picard_processed/3.1_sortSam/ -p
mkdir $analysis_dir/3_Picard_processed/3.2_AddOrReplaceReadGroups -p
mkdir $analysis_dir/3_Picard_processed/3.3_markDuplicates -p
mkdir $analysis_dir/3_Picard_processed/3.3_markDuplicates/deduppedBamSortedFiles -p
mkdir $analysis_dir/3_Picard_processed/temp -p
#mkdir $analysis_dir/sam_files/ -p

# After mapping, the aligned sam files are scattered in sample-specific folders
# We want to find all .sam files and move them to one directory
# find "$sam_files" -maxdepth 2 -type f -name  "*.sam"  -printf "/%P\n" | while read FILE ; do filename=$(basename $FILE) ; mv "$sam_files/$FILE" "$analysis_dir/sam_files/$filename"; done 
# files=$(echo $analysis_dir/sam_files/*.sam | xargs -n 1 basename)

# 1. Sorting
#        --SORT_ORDER: coordinate or queryname (required)
#       at https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-



# 2. AddOrReplaceReadGroups
#       --RGLB: Read-Group library (required)
#       --RGPL: Read-Group platform, e.g. ILLUMINA, SOLID. (required)
#       --RGPU: Read-Group platform unit. (required)
#       --RGSM: Read-Group sample name (required)
#       https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-


#        --MAX_RECORDS_IN_RAM: 500000 (default value; optional): When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed
#       --COMPRESSION_LEVEL: Compression level for all compressed files created (e.g. BAM and VCF). Default: 5. Can set to 0 (no compression) ---> faster, as no time is spent to compress the file. But takes more disk space.
#       --VALIDATION_STRINGENCY: STRICT or SILENT. Validation stringency for all SAM files read by this program. 


     java  $memory -Djava.io.tmpdir=$analysis_dir/3_Picard_processed/temp -jar $sw/picard.jar SortSam -I $sam_files/$samplename/$filename -O $analysis_dir/3_Picard_processed/3.1_sortSam/$samplename.ordered.sam -SORT_ORDER coordinate -TMP_DIR $analysis_dir/3_Picard_processed/temp -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/3_Picard_processed/3.1_sortSam/3.1.picard_sorting.$samplename.txt 2>&1
    
     java  $memory -Djava.io.tmpdir=$analysis_dir/3_Picard_processed/temp -jar $sw/picard.jar AddOrReplaceReadGroups -I $analysis_dir/3_Picard_processed/3.1_sortSam/$samplename.ordered.sam -O $analysis_dir/3_Picard_processed/3.2_AddOrReplaceReadGroups/$samplename.addReplace.bam -RGLB $filename -RGPL ILLUMINA -RGPU ILLUMINA -RGSM $filename -RGID $samplename -SORT_ORDER coordinate -COMPRESSION_LEVEL 5 -VALIDATION_STRINGENCY SILENT -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/3_Picard_processed/3.2_AddOrReplaceReadGroups/3.2.AddOrReplaceReadGroups.$samplename.txt 2>&1
     
     java  $memory -Djava.io.tmpdir=$analysis_dir/3_Picard_processed/temp -jar  $sw/picard.jar MarkDuplicates -I $analysis_dir/3_Picard_processed/3.2_AddOrReplaceReadGroups/$samplename.addReplace.bam -O $analysis_dir/3_Picard_processed/3.3_markDuplicates/$samplename.marked.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -ASSUME_SORT_ORDER coordinate -M $analysis_dir/3_Picard_processed/3.3_markDuplicates/deduppedBamSortedFiles/$samplename.output.metrics -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/3_Picard_processed/3.3_markDuplicates/3.3_markDuplicates.$samplename.txt 2>&1
     

