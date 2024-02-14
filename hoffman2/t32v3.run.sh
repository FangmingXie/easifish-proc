#!/bin/bash

# use current wd
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=20G
#$ -pe shared 2 
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3 
conda activate napari

python3 ./32.correction_s0_v3.py