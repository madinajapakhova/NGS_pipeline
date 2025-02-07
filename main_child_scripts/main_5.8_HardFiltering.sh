#!/bin/bash
set -e
analysis_dir=.../analysis_folder/...

echo "Started at `date`"

        sbatch -J 5.8_HardFiltering ./main_child_scripts/child_5.8_HardFiltering.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


