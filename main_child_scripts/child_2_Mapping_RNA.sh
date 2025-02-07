#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=2_child_Mapping # write a meaningful name
#SBATCH --time=2-00:00:00 # dd-hh:mm:ss
#SBATCH --output=.../logs/%x-%j.log 
#SBATCH --error=.../%x-%j.err 
#SBATCH --mem-per-cpu=300000  # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user example@whatever.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type=FAIL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge # to prevent conflicts 
module load apps/star/2.7.9a

analysis_dir=$1 
processed_files=$2                                      
paired_end=$3 
filename_extention=$4 
samplename=$5 
readname=$6 
genomeDir=$7
reference_genome38=$8
reference_gtf=$9
threads=${10}
overhang=${11}

echo 'Started'

./2_Mapping_RNA.sh $analysis_dir $processed_files $paired_end $filename_extention $samplename $readname $genomeDir $reference_genome38 $reference_gtf $threads $overhang

echo "Finished at `date`"
echo "__DONE__"
