
#bash template
#!/bin/bash
#SBATCH --time=3-00:00:00   # Time limit for the job (REQUIRED).
#SBATCH --job-name=ra2263   # Job name
##SBATCH --gres=gpu:1      #number of gpu per task
#SBATCH --nodes=1           # number of nodes
##SBATCH --ntasks=1          # Number of cores for the job. Same as SBATCH -n 1
##SBATCH --mem=32g           #32 GB ram asked
##SBATCH --partition=CAL48M192_L #choose one of the nodes
##SBATCH --partition=CAC48M192_L #choose one of the nodes
#SBATCH --partition=SKY32M192_L
##SBATCH --partition=P4V16_HAS16M128_L
##SBATCH --partition=P4V12_SKY32M192_L
#SBATCH -e ra2265-%j.err
#SBATCH -o ra2265-%j.log
#SBATCH -A col_rma422_uksr  # Project allocation account name (REQUIRED)
#SBATCH --mail-type FAIL    # Send email when job starts/ends
#SBATCH --mail-user tpbh222@uky.edu# Where email is sent to (optional)


echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "
echo "`date`"

#loading the essential modules
module purge
module load cmake/3.14.3
module load ccs/conda/root-6.22.0+geant4-10.6.0 
$HOME/gesim/gesim ra2265.mac          #gesim
echo "DONE!"
echo "`date`"
