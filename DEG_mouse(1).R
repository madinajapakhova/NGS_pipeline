#!/usr/bin/Rscript 
# srun --pty -c 20 --mem 70GB --time 12:00:00 bash

# R version used for this script: 4.1.1 from singularity
# To launch it please run the following two lines:
# module load apps/singularity/3.9.2
# singularity exec /group/bioinf/Software/singularity/R441.sif R

# ******************************************************************************
# Load packages 
# ******************************************************************************
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(DESeq2)
library(ape)
library(biomaRt) 
library(AcidGenomes) # for creating a GRanges object on a gene level
library(txdbmaker) # for function makeGRangesFromGff from AcidGenomes
library(ggplot2)
library(EnhancedVolcano)
library(stringr) # needed this to change column names in a count matrix based on a pattern
library(dplyr) # to manage dataframes
library(readxl)
library(gplots)
library(scatterplot3d) 
library(pheatmap)
library(preprocessCore)   
# ******************************************************************************

rm(list=ls())
options(timeout = 120000) # To download a large file (> 50 MB) that it exceeds the "timeout" option built-into R; miliseconds
working_dir <- "/group/bioinf_biomarkers_rna/madina/mouse_check"
setwd(working_dir)

# ******************************************************************************
# PLEASE PROVIDE: 1) SAMPLE DATA (CLINICAL DATA/METADATA) 
#                 2) PATH TO THE FOLDER WHERE YOU KEEP CURATED BAM FILES READY FOR ANALYSIS
#                 3) GTF FILE. USE THE SAME GTF FILE THAT WAS USED FOR MAPPING
# ******************************************************************************
metadata <- read_excel("phenotype_data_STK11.xlsx") # samplenames - column 1

# Read quantification takes some time. The more samples - the more time. So if you already have the count matrix, just load it as a .csv file
count_matrix_generated <- "No" # set "Yes" if you already have the count matrix

if (count_matrix_generated == "Yes") {
    read_counts_0 <- read.csv("/group/bioinf_biomarkers_rna/mice_stk11_christian/Ballal/Generated_data/read_counts_raw.csv", header=TRUE, row.names=1)
} else {
    path_to_bam_files <- "/group/bioinf_biomarkers_rna/mice_stk11_christian/Ballal/STK11_results/3_Picard_processed/3.3_markDuplicates"
    #gtfFile <- "/group/bioinf_biomarkers_rna/mice_stk11_christian/Ballal/reference_genome/GRCm39_GTF_GFF/genomic.gtf"
    gtfFile <- "/group/bioinf_biomarkers_rna/mice_stk11_christian/Ballal/reference_genome/new/mm39.ncbiRefSeq.gtf"
}
# ******************************************************************************

# ******************************************************************************
# PHENOTYPE DEFINITION
# ******************************************************************************
reference <- "positive"
metadata <- subset(metadata, (phenotype == "negative" | phenotype == "positive"))
metadata$phenotype <- factor(metadata$phenotype)
metadata <- metadata[!is.na(metadata$phenotype), ] 

# ******************************************************************************
# Initialize the backend for parallelising. (BiocParallel package)
# Defaults to all cores available as determined by detectCores
# ******************************************************************************
use_all_cores <- "Yes"
ncores <- 12 # number of workers

if (use_all_cores == "Yes") {
    multicoreParam <- MulticoreParam()
    register(multicoreParam)
    registered()
} else {
    multicoreParam <- MulticoreParam(ncores=ncores)
    register(multicoreParam)
    registered()
}

# ******************************************************************************
# A FEW WORDS ON THE COUNT MATRIX
# ******************************************************************************
# DeSeq2 requires as input an unnormalized count matrix, therefore we first need to generate a count matrix
# Count matrix: how many reads overlap with a genomic region (feature, e.g. gene; transcript; exon) in each sample 
# Why unnormalized? DeSeq2 internally corrects for library size

# Removing cigar strings from bam files is not needed

# There are multiple ways to generate a count matrix
# We will use SummarizeOverlaps from Bioconductor
# ******************************************************************************

# ******************************************************************************
# GENERATING THE COUNT MATRIX
# ******************************************************************************

# If count matrix has previously been generated, skip this quantification step directly go to checking COMPATIBILITY OF COUNT MATRIX AND SAMPLE DATA 

if (count_matrix_generated == "Yes") {
    print("Count matrix has already been generated. Go to checking COMPATIBILITY OF COUNT MATRIX AND SAMPLE DATA.")
} else {
    # SummarizeOverlaps requires as input two things: bam files (generated with our pipeline) and features (genomic regions)
    # Input 1 for SummarizeOverlaps: bam files 
    filenames <- dir(path_to_bam_files,pattern=".bam$",full.names=TRUE)
    bfl <- BamFileList(filenames, yieldSize = 50000, index = character()) 

    # Input 2 for SummarizeOverlaps: features
    # Features can be genes; transcripts; exons ---> discussion
    # Features are stored as a GRanges or GRangesList object
    # GRanges: a genomic region is described using one range, e.g. One region=one gene
    # GRangesList: a genomic region is described using several ranges, e.g. One region={gene; transcript}
    # In this script we do gene expression analysis on a gene level 

    # Creating a GRanges object on a gene level using a gene model provided by the GTF file 
    # Tool: makeGRangesFromGFF from AcidGenomes library
    #genedb<-makeGRangesFromGFF(gtfFile, level="genes") # Warning message: reported and downloaded sizes are not same (?) Probably caused by that the file size exceeds R's timeout option
    
    genedb<-makeGRangesFromGff(gtfFile,level = c("genes"), ignoreVersion = FALSE, extraMcols = TRUE)


    # Generating count matrix
    counteByg <- bplapply(
        bfl, 
        function(x) summarizeOverlaps(
        genedb,
        x,
        mode = "Union", 
        ignore.strand = TRUE, 
        inter.feature = TRUE,
        singleEnd = TRUE, 
        BPPARAM = multicoreParam))
    
    countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
    rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
    colnames(countDFeByg) <- names(bfl)
    read_counts_0 <- as.data.frame(countDFeByg)
    names(read_counts_0) <- substr(names(read_counts_0), start=1, stop=6) # Keep only first 6 chr. in sample name for future simplicity
}

# Count matrix on a gene level has been generated 

# 

if (count_matrix_generated == "Yes") {
    print("Count matrix has already been generated. Make sure you have checked COMPATIBILITY OF COUNT MATRIX AND SAMPLE DATA.")
} else {
    write.csv(read_counts_0, file="read_counts_raw_stk11.csv")
}

# ******************************************************************************
# COMPATIBILITY OF METADATA AND READ COUNTS DATAFRAME
# ******************************************************************************
metadata$STK11_Data[which(metadata$STK11_Data == "248644-TB-WT")] <- "248644"
metadata$STK11_Data[which(metadata$STK11_Data == "248645-TB-NC1.4")] <- "248645"
metadata$STK11_Data[which(metadata$STK11_Data == "248646-TB-NC2.2")] <- "248646"
metadata$STK11_Data[which(metadata$STK11_Data == "248649-TB-Stk11-KO1.5")] <- "248649"
metadata$STK11_Data[which(metadata$STK11_Data == "248650-TB-Stk11-KO2.2")] <- "248650"
metadata$STK11_Data[which(metadata$STK11_Data == "261861-TB-Stk11-KO1.8-N-trimmed")] <- "261861"

rownames(metadata) <- metadata$STK11_Data

# Check if all of the rownames of metadata are contained as colnames in count data
compatibility_check_1 <- all(rownames(metadata) %in% colnames(read_counts_0)) 
# Check if rownames of metadata are same as colnames as count data
compatibility_check_2 <- all(rownames(metadata) == colnames(read_counts_0)) 
# Check if first column of metadata corresponds to colnames of count data
compatibility_check_3 <- all(metadata[,1]==colnames(read_counts_0))

if (compatibility_check_1 == TRUE & compatibility_check_2 == TRUE & compatibility_check_3 == TRUE) {
    print("Count matrix and sample data are compatible. All good.")
    read_counts <- read_counts[, unique(rownames(metadata))]
} else {
    stop("Count matrix and sample data are not compatible. Please make them compatible")
}

read_counts_0 <- read_counts_0[, unique(rownames(metadata))]

# Check if all of the rownames of metadata are contained as colnames in count data
compatibility_check_1 <- all(rownames(metadata) %in% colnames(read_counts_0)) 
# Check if rownames of metadata are same as colnames as count data
compatibility_check_2 <- all(rownames(metadata) == colnames(read_counts_0)) 
# Check if first column of metadata corresponds to colnames of count data
compatibility_check_3 <- all(metadata[,1]==colnames(read_counts_0))

if (compatibility_check_1 == TRUE & compatibility_check_2 == TRUE & compatibility_check_3 == TRUE) {
    print("Count matrix and sample data are compatible. All good.")
    read_counts_0 <- read_counts_0[, unique(rownames(metadata))]
} else {
    stop("Count matrix and sample data are not compatible. Please make them compatible")
}

# ******************************************************************************
# Differential gene expression analysis
# ******************************************************************************
dds <- DESeqDataSetFromMatrix(countData = read_counts_0,
                              colData = metadata,
                              design = ~ phenotype)

# Pre-filtering genes with low counts (optional). For visualisations and speeding up processing. 
# Filter out genes where there are less than 1 sample with counts greater than or equal to 0.
smallestGroupSize <- 1 # number of samples
keep <- rowSums(counts(dds) >= 0) >= smallestGroupSize
dds <- dds[keep,]

# Note: by default the control level will be the first level or the one that comes first alphabetically. To control which level is control and which one is treatment we use function relevel later. 
# specifying which group is the reference: relevel 
dds$phenotype <- relevel(dds$phenotype, ref = reference)

# Check for size factors     
# Comparing each sample to an "average" pseudosample through a rayio. DeSeq2 uses geometric mean, because it is more stable to outliers than arithmetic mean
# Should be around 1.                     
dds <- estimateSizeFactors(dds)  
sizeFactors(dds) # Good values are around 1, meaning most samples are similar to average pseudosample.
#write.csv(as.data.frame(sizeFactors(dds)), file="SF_estimates_pdac_cp.csv")
                          
# run differential expression analysis  
# DESeq, results, and lfcShrink can be parallelized                    
ddsRNASeq <- DESeq(dds, parallel=TRUE, BPPARAM = multicoreParam)

# extract DESeq normalized counts for PCA
deseq_counts <- counts(dds, normalized=T)
#vsd_counts <- vst(dds, blind=FALSE)
#write.csv(deseq_counts, file="deseq_normalized_counts.csv")

ddsRNASeq <- results(ddsRNASeq, parallel=TRUE, BPPARAM = multicoreParam)
summary(ddsRNASeq)

# Sort this table by p-value (smaller p-values on top)
orderedRes <- ddsRNASeq[ order(ddsRNASeq$padj), ]
orderedRes <- orderedRes[!(orderedRes$baseMean == 0 & is.na(orderedRes$log2FoldChange)), ]
write.csv(as.data.frame(orderedRes), file="DESeq_results.csv") # annotate this table






































































