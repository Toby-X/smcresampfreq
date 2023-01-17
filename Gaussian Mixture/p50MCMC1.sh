#!/bin/bash
#SBATCH -J SMCS_5001
#SBATCH -n 2
#SBATCH -t 72:00:00
#SBATCH -N 2
#SBATCH -c 48
#SBATCH -p queue
#SBATCH -o /public1/home/scf0347/ResampFreq/GaussianMixture/p50MCMC1.txt
source /public1/soft/modules/module.sh
module load anaconda
source activate R4.2.2
Rscript /public1/home/scf0347/ResampFreq/GaussianMixture/p50_MCMC_1_par.R