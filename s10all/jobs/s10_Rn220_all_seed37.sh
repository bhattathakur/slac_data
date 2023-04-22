#!/bin/bash
#SBATCH --account=shared --partition=shared
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5g
#SBATCH --time=1-10:10:00
#SBATCH --error %j.err
#SBATCH --output %j.log
#SBATCH --mail-type FAIL                        # Send email when job starts/ends
#SBATCH --mail-user tpbh222@uky.edu            # Where email is sent to (optional)


singularity exec --bind /scratch/t/tbhatta --home /sdf/home/t/tbhatta /sdf/group/nexo/nexo_spack18_g4.7.2_latest.sif ./s10_Rn220_all_seed37submit.sh
