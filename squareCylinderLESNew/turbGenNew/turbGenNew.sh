#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N TURB_GEN
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe mpi 72
#$ -P quanah

module load matlab/R2017b
matlab -nodisplay -nosplash < main.m
#matlab -nodisplay -nosplash < post_inhomo.m

