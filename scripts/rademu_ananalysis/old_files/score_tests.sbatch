#!/bin/bash
#SBATCH --job-name hdw   
#SBATCH --partition short        
#SBATCH --ntasks 1                   
#SBATCH --time 0-11:50                
#SBATCH --mem-per-cpu=800       
#SBATCH --array=1-272                 
#SBATCH -o hdw_%j.out           
#SBATCH -e hdw_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript hdw/cluster_score_test_script.R $SLURM_ARRAY_TASK_ID