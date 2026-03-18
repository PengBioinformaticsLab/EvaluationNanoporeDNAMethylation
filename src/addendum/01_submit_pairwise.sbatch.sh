#!/bin/bash
#SBATCH --job-name=pairwise_EM_ONT
#SBATCH --account=r00302
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stevbroo@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

set -euo pipefail

echo "Running job: $SLURM_JOB_NAME"

module load r

Rscript --vanilla addendum/pairwise_pearsons_EM_ONT.R
Rscript --vanilla addendum/pairwise_pearsons_ONT_blood.R
Rscript --vanilla addendum/pairwise_pearsons_EPIC_blood.R