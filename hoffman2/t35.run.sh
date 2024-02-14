#!/bin/bash

# use current wd
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=23:59:00,h_data=20G
#$ -pe shared 2 
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3 
conda activate napari

python3 ./35.correction_s1_lt185_r2_batch2.py