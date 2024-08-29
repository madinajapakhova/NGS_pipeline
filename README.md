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
[**GenomicAlignments**](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html) - for working with aligned samples, specifically - counting reads     
[**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) - for differential gene expression analysis   
[**AcidGenomes**](https://github.com/acidgenomics/r-acidgenomes) - for working with annotated genomes        
[**biomaRt**](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) - for retrieving gene annotations from different databases      
[**BiocParallel**](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) - for parallelizing some of Bioconductor functionalities           
[**EnhancedVolcano**](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) - for plotting results of gene expression analysis           

### Other R packages      
[**stringr**] - for working with strings     
[**dplyr**] - for managing dataframes       

### Singularity container    
We will provide all tools exploited by the pipeline, including $R$ with all necessary packages, installed within a singularity container. Therefore, you can download the container and execute the pipeline with all functions within the container. You only need to have [**Singularity**](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) installed.  

## Pipeline processes    
![](https://github.com/madinajapakhova/NGS_pipeline/blob/main/pipeline_processes.png)   

## Parallelization and execution   

### Main-child-pipeline design    
All pipeline processes, except FastQC*, are executed independently per each sample in a parallel fashion on HPC. Parallelization is achieved via means of **main-child-pipeline** 3-script design. The **main** script runs the processe per a list of samples, the **child** script contains job description, i.e. how much memory to allocate, and the **pipeline** script is a script with one of the pipeline processes. **You only need to modify the *main* script** by providing the expected input, e.g. path to files.  The **main** script then calls the **child**, and the **child** then calls the **pipeline** process script (e.g. data curation). The **child** and **pipeline** scripts have to stay unaffected!   

* Scripts using FastQC have to be used on their own - without the **main** and **child** scripts. The reason is that FastQC has a nice in-built support for independent parallel processing. It is simple and efficient, therefore the pipeline relies on the FastQC in-built parallelization. 

### HPC   
The pipeline is executed on Linux-based HPC via the [**Slurm workload manager**](https://slurm.schedmd.com/sbatch.html) to allocate jobs. Make sure your computing facility supports [Slurm](https://slurm.schedmd.com/sbatch.html).             

