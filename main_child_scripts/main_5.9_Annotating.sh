#!/bin/bash
set -e
analysis_dir=.../analysis_folder/...

echo "Started at `date`"

        sbatch -J 5.9_Annotating ./main_child_scripts/child_5.9_Annotating.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


