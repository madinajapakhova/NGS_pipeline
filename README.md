# NGS pipeline: [the repository is currently under development]
We present our pipeline for processing NGS data, optimized for parallel execution on HPC systems. OS compatibility: Linux. 

![](https://github.com/madinajapakhova/NGS_pipeline/blob/main/workflow_overview.png)  

## Tools 
### Bioinformatics software    
[**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - a Java-based software that generates html reports evaluating the quality of FASTQ data.        
[**MultiQC**](https://multiqc.info/) - aggregates analysis reports from multiple samples into a unified report       
[**UMI-tools**](https://umi-tools.readthedocs.io/en/latest/) - handles Unique Molecular Identifiers (UMIs)               
[**Cutadapt**](https://cutadapt.readthedocs.io/en/stable/) - removes unwated sequences (e.g. adapters, primers, poly-A tails) from the sequenced reads            
[**Trimmomatic**](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) - trims sequences of low quality, as well as adapters, primers and other unwanted sequences from Illumina (FASTQ) data             
[**STAR**](https://github.com/alexdobin/STAR) - aligns RNA-seq to the reference genome       
[**BWA**](https://bio-bwa.sourceforge.net/) - aligns DNA sequences to the reference genome      
[**Qualimap**](http://qualimap.conesalab.org/) - evaluates the quality of aligned samples        
[**GATK**](https://gatk.broadinstitute.org/hc/en-us) - performs variant calling     

### Bioconductor R packages
[**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) - Bioconductor R package for differential gene expression analysis       
[To be filled]

## Pipeline processes    
![](https://github.com/madinajapakhova/NGS_pipeline/blob/main/pipeline_processes.png) 
