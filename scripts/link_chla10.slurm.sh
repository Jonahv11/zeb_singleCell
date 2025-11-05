#! /bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-ewing_beti-sponsored
#SBATCH --mem-per-cpu=40g
#SBATCH --time=0-4:00:00
#SBATCH --job-name=linkpeaks_chla10
#SBATCH --chdir=/data/hps/assoc/private/ewing_beti/user/jvale8/slurm_logs/zeb_singleCell/links/
#SBATCH --mail-user=jonah.valenti@seattlechildrens.org
#SBATCH --mail-type=ALL

# to run script: sbatch <scriptname>.slurm.sh

cd "/data/hps/assoc/private/ewing_beti/user/jvale8/zeb_singleCell/scripts/"

apptainer exec --bind /data/hps/assoc /data/hps/assoc/public/bioinformatics/container/posit/posit-base-20250818.sif Rscript link_chla10.R

