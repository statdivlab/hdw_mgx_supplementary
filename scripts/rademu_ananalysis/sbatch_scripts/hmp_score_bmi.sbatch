#!/bin/bash
#SBATCH --job-name hmp_bmi   
#SBATCH --partition medium       
#SBATCH --ntasks 1                   
#SBATCH --time 3-11:50             
#SBATCH --mem-per-cpu=1200       
#SBATCH --array=1-272                 
#SBATCH -o hmp_score_bmi_%j.out           
#SBATCH -e hmp_score_bmi_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript hdw/scripts/hmp_score_tests_cluster_bmi.R $SLURM_ARRAY_TASK_ID