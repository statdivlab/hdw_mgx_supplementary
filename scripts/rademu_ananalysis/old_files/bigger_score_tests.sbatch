#!/bin/bash
#SBATCH --job-name hdw_bigger   
#SBATCH --partition medium       
#SBATCH --ntasks 1                   
#SBATCH --time 6-11:50                
#SBATCH --mem-per-cpu=1600       
#SBATCH --array=1-412                 
#SBATCH -o hdw_bigger_%j.out           
#SBATCH -e hdw_bigger_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript hdw/bigger_score_tests.R $SLURM_ARRAY_TASK_ID