#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=5.4_5.6_Genotyping # write a meaningful name
#SBATCH --time=12:00:00 # dd-hh:mm:ss
#SBATCH --output=.../logs/%x-%j.log # just write your zih name instead of sxxxxx, and keep the rest as they are 
#SBATCH --error=.../logs/%x-%j.err  # just write your zih name instead of sxxxxx,  and keep the rest as they are 
#SBATCH --mem-per-cpu=200000 # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user example@whatever.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type=FAIL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge 
module load apps/java/8u41-SE
module load apps/samtools/1.12
module load apps/singularity/3.7.3

analysis_dir=$1 
chrom=$2 

echo 'Started'

./5.4_5.6_Genotyping.sh $analysis_dir $chrom 


echo "Finished at `date`"
echo "__DONE__"
