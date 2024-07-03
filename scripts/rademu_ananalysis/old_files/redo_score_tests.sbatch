#!/bin/bash
#SBATCH --job-name hdw_redo   
#SBATCH --partition medium      
#SBATCH --ntasks 1                   
#SBATCH --time 3-11:50                
#SBATCH --mem-per-cpu=1600       
#SBATCH --array=1-27                 
#SBATCH -o hdw_redo_%j.out           
#SBATCH -e hdw_redo_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript hdw/redo_score_test.R $SLURM_ARRAY_TASK_ID