#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=2.1_GeneratingRefGenome # write a meaningful name
#SBATCH --time=06:00:00 # days-hh:mm:ss
#SBATCH --output=.../logs/%x-%j.log 
#SBATCH --error=.../%x-%j.err  
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
module purge 
module load apps/star/2.7.9a

analysis_dir=$1 
threads=$2
reference_genome38=$3
reference_gtf=$4 
temp_folder=$5 
                           
echo 'Started'

./2.1_GeneratingRefGenome.sh $analysis_dir $threads $reference_genome38 $reference_gtf $temp_folder 



echo "Finished at `date`"
echo "__DONE__"
