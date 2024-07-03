#!/bin/bash
#SBATCH --job-name hmp   
#SBATCH --partition medium       
#SBATCH --ntasks 1                   
#SBATCH --time 3-11:50             
#SBATCH --mem-per-cpu=800       
#SBATCH --array=1-272                 
#SBATCH -o hmp_score_%j.out           
#SBATCH -e hmp_score_%j.err            
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=teichs@uw.edu 
module load R
Rscript hdw/scripts/hmp_score_tests_cluster.R $SLURM_ARRAY_TASK_ID