#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=5.1_HaplotypeCaller.sh # write a meaningful name
#SBATCH --time=24:00:00 # 100 hours as an example
#SBATCH --output=/home/bioinf/maja488d/logs/%x-%j.log # just write your zih name instead of sxxxxx, and keep the rest as they are 
#SBATCH --error=/home/bioinf/maja488d/logs/%x-%j.err  # just write your zih name instead of sxxxxx,  and keep the rest as they are 
#SBATCH --mem-per-cpu=100000 # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user madina.japakhova@mailbox.tu-dresden.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type=FAIL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge # to prevent conflicts 
module load apps/java/8u41-SE

analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1 main analysis folder where all performed  steps of the pipeline are saved
bam_files=$2 # /group/bioinf/Data/Madina_test/pipeline/v1/Picard_processed/3.3_markDuplicates/                                          
filename_extention=$3 # .bam
samplename=$4 # name of the sample, e.g. DNA1
echo 'Started'

./E.Coli_variant_calling/5.1._HaplotypeCaller.sh $analysis_dir $bam_files $filename_extention $samplename 



echo "Finished at `date`"
echo "__DONE__"
