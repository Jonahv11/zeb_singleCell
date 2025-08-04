#! /bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --partition=cpu-core-sponsored
#SBATCH --account=cpu-ewing_beti-sponsored
#SBATCH --mem-per-cpu=60g
#SBATCH --time=0-3:00:00
#SBATCH --job-name=tgfb1_treated_multiome_data_processing
#SBATCH --chdir=/data/hps/assoc/private/ewing_beti/user/jvale8/slurm_logs/zeb_singleCell/tgfb1_treated/
#SBATCH --mail-user=jonah.valenti@seattlechildrens.org
#SBATCH --mail-type=ALL

cd "/data/hps/assoc/private/ewing_beti/user/jvale8/zeb_singleCell/cell_line_contrasts/src/"

apptainer exec --bind /data/hps/assoc /data/hps/assoc/public/bioinformatics/container/posit/posit-base-20250605.sif Rscript tgfb1_processing.R