#!/bin/bash
#SBATCH --partition=batch
#SBATCH --qos=240c-1h_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name='TrajAnalysis'
#SBATCH --output=TrajAna.%J.out
#SBATCH --error=TrajAna.%J.err
#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
#SBATCH --mail-type=ALL
module load gromacs/2021

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# MAIN
WD=$SLURM_SUBMIT_DIR
cd $WD
# PREPARATION PART WAS REDACTED 

#perform RMSF calculations for ligand
echo "Performing Ligand RMSF Calculations..."
if ! [ -f *RMSF* ]; then
	echo "21" | gmx rmsf -f xxxx_mdout_center.xtc -s xxxx_minimization.tpr -o xxxx_RMSF.xvg -n index.ndx
else
	echo "RMSF file has been found. Skipping this step..."
fi

#perform gyration analysis

#Calculating Complex Gyration
if ! [ -f *ProtLig_gyrate* ]; then
	echo "19" | gmx gyrate -f xxxx_mdout_center.xtc -s xxxx_minimization.tpr -n index.ndx -o xxxx_ProtLig_gyrate.xvg
else
	echo "Gyration Data has been found. Skipping this step..."
fi

#Calculating Protein Gyration
if ! [ -f *ProtOnly_gyrate* ]; then
	echo "1" | gmx gyrate -f xxxx_mdout_center.xtc -s xxxx_minimization.tpr -n index.ndx -o xxxx_ProtOnly_gyrate.xvg
else
	echo "Gyration Data has been found. Skipping this step..."
fi

#get the hydrogen bond distribution
if [[ ! -f hbnum.xvg ]] && [[ ! -f hbdist.xvg ]]; then
	echo "1 13" | gmx hbond -f xxxx_mdout_center.xtc -s xxxx_minimization.tpr -n index.ndx -num -dist
else
	echo "Hydrogen Bond Distribution data has been found. Skipping this step..."
fi
