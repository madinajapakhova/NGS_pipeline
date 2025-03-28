#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=child_3_Picard.sh # write a meaningful name
#SBATCH --time=12:00:00 # dd-hh:mm:ss
#SBATCH --output=.../logs/%x-%j.log # just write your zih name instead of sxxxxx, and keep the rest as they are 
#SBATCH --error=.../logs/%x-%j.err  # just write your zih name instead of sxxxxx,  and keep the rest as they are 
#SBATCH --mem-per-cpu=200000  # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user example@whatever.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type ALL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge # to prevent conflicts 
module load apps/singularity/3.7.3
module load apps/java/20.0.1

analysis_dir=$1 
sam_files=$2
filename_extention=$3 
samplename=$4 
                                           
echo 'Started'

./3_PicardProcessing.sh $analysis_dir $sam_files $filename_extention $samplename  



echo "Finished at `date`"
echo "__DONE__"
