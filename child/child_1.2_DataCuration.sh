#!/bin/bash
#SBATCH -A pipeline
#SBATCH --job-name=child_1.2_DataCuration.sh # write a meaningful name
#SBATCH --time=1-12:00:00 # dd-hh:mm:ss "1-12:00:00" means "one day and twelve hours"
#SBATCH --output=.../logs/%x-%j.log 
#SBATCH --error=.../%x-%j.err  
#SBATCH --mem-per-cpu=5000  # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=1 # number of the request cpu if you have parallel processing
#SBATCH --mail-user example@whatever.de	### tell the batch system your email address to get updates about the jobs status
#SBATCH --mail-type=FAIL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 

# to see the available modules type this in terminal: module avail

## modules
module purge 
module load apps/cutadapt/3.4
module load apps/java/16.0.1
module load apps/python3/3.9.5

analysis_dir=$1 
folder_nextstep=$2                                         
paired_end=$3 
filename_list_exclude=$4   
softwares=$5 
UMI_extraction=$6 
UMI_pattern=$7 
UMI_pattern2=$8 
umi_method=$9
cut_adapters=${10} 
trimming_software=${11} 
adapter_fwd=${12} 
adapter_rev=${13} 
illumina_adapters_file_PE=${14} 
illumina_adapters_file_SE=${15} 
low_quality=${16} 
remove_poly=${17}
threads=${18} 
samplename=${19}
filename=${20}
filename2=${21}

echo 'Started'

./scripts/1.2_DataCuration.sh $analysis_dir $folder_nextstep $paired_end $filename_list_exclude $softwares $UMI_extraction $UMI_pattern $UMI_pattern2 $umi_method $cut_adapters $trimming_software $adapter_fwd $adapter_rev $illumina_adapters_file_PE $illumina_adapters_file_SE $low_quality $remove_poly $threads $samplename $filename $filename2



echo "Finished at `date`"
echo "__DONE__"
