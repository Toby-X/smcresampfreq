#!/bin/bash
#SBATCH -J Res_10001
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -p queue
#SBATCH -o /public1/home/scf0347/ResampFreq/GaussianMixture/ResResample/p100MCMC1.txt
source /public1/soft/modules/module.sh
module load anaconda
source activate R4.2.2
Rscript /public1/home/scf0347/ResampFreq/GaussianMixture/ResResample/p100_MCMC1_res.R