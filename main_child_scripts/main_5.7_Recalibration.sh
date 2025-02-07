#!/bin/bash
set -e
analysis_dir=.../analysis_dir/...

echo "Started at `date`"

        sbatch -J 5.7_Recalibration ./main_child_scripts/child_5.7_VQSR.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


