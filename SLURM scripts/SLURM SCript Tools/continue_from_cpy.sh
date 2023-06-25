#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --qos=12c-1h_2gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --exclude=saliksik-gpu-02
#SBATCH --job-name='com_Pr1_rrun'
#SBATCH --output=JobName.%J.out
#SBATCH --error=JobName.%J.err
#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
#SBATCH --mail-type=ALL
module load gromacs/2021.4_cuda-11.4.3


echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# Set stack size to unlimited
ulimit -s unlimited

# MAIN
WD=$SLURM_SUBMIT_DIR
cd $WD


#takenote that -deffnm are the same with last run
# completion at 50000000 



if [ -f "$WD/com_Pr1_E484A_mdout.log" ]; then
	if grep -Fq "Writing checkpoint, step 50000000" $WD/com_Pr1_E484A_mdout.log; then
		exit
	else
		gmx mdrun -deffnm com_Pr1_E484A_mdout -s $WD/com_Pr1_E484A_mdout.tpr -x $WD/com_Pr1_E484A_mdout.xtc  -cpi $WD/com_Pr1_E484A_mdout.cpt -maxh 72
		sbatch continue_from_cpy.sh
		exit
	fi
fi
	
