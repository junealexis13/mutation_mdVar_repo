#!/bin/bash

WSUBMIT=$1
WORKCODE=$2

if [[ $WSUBMIT == "gpu" ]]; then
	#SBATCH --partition=gpu
	#SBATCH --qos=12c-1h_2gpu
	#SBATCH --nodes=2
	#SBATCH --ntasks-per-node=8
	#SBATCH --exclude=saliksik-gpu-02
	#SBATCH --job-name="$WORKCODE"
	#SBATCH --output=JobName.%J.out
	#SBATCH --error=JobName.%J.err
	#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
	#SBATCH --mail-type=ALL

	module load anaconda
	mamba activate gromacs-2021.3-cuda



elif [[ $WSUBMIT == "cpu" ]]; then
	#SBATCH --partition=batch
	#SBATCH --qos=240c-1h_batch
	#SBATCH --nodes=2
	#SBATCH --ntasks-per-node=32
	#SBATCH --job-name="$WORKCODE"
	#SBATCH --output=JobName.%J.out
	#SBATCH --error=JobName.%J.err
	#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
	#SBATCH --mail-type=ALL

	module load anaconda
	mamba activate gromacs-2021.3-cuda



else
	echo "Mode not specified"
	exit


echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# MAIN
WD=$SLURM_SUBMIT_DIR
cd $WD
# PREPARATION PART WAS REDACTED 


if [[ -f $WD/nvt_omic.mdp ]] && [[ -f $WD/xxxx_min.gro ]]; then
        #NVT equilibration
        gmx grompp -f $WD/nvt_omic.mdp -c $WD/xxxx_min.gro -r $WD/xxxx_min.gro -p $WD/topol.top -n $WD/xxxx_min.ndx -o $WD/xxxx_equi.tpr
        gmx mdrun -deffnm xxxx_nvt -s $WD/xxxx_equi.tpr -x $WD/xxxx_nvt.xtc
fi

if [[ -f $WD/npt_omic.mdp ]] && [[ -f $WD/xxxx_nvt.gro ]]; then
        #NPT Equilibration
        gmx grompp -f $WD/npt_omic.mdp -c $WD/xxxx_nvt.gro -r $WD/xxxx_nvt.gro -t $WD/xxxx_nvt.cpt -o $WD/xxxx_npt.tpr -p $WD/topol.top -maxwarn 3
        gmx mdrun -deffnm xxxx_npt -s $WD/xxxx_npt.tpr -x $WD/xxxx_npt.xtc

        #Getting NPT            energy log
        #echo "18 0" | gmx energy -f $WD/complex_lig_prot_npt.edr -o $WD/complex_lig_prot_pressure.xvg

        #Getting NPT density log
        #echo "24 0" | gmx energy -f $WD/complex_lig_prot_npt.edr -o $WD/density.xvg

fi

if [[ -f $WD/trial_md.mdp ]] && [[ -f xxxx_npt.gro ]]; then
        #Executing MD RUN
	if [[ $WSUBMIT == "cpu" ]]; then
        	gmx grompp -f $WD/trial_md.mdp -c $WD/xxxx_npt.gro -t $WD/xxxx_npt.cpt -p $WD/topol.top -o $WD/xxxx_mdout.tpr
        	gmx mdrun -ntmpi 16 -ntomp 5 -npme 4 -deffnm xxxx_mdout -s $WD/xxxx_mdout.tpr -x $WD/xxxx_mdout.xtc
	elif [[ $WSUBMIT == "gpu" ]]; then
		gmx grompp -f $WD/trial_md.mdp -c $WD/xxxx_npt.gro -t $WD/xxxx_npt.cpt -p $WD/topol.top -o $WD/xxxx_mdout.tpr
		gmx mdrun -deffnm production -nt 6 -ntmpi 2 -ntomp 3 -npme 1 -nb gpu -pme cpu -bonded gpu -s $WD/xxxx_mdout.tpr -x $WD/xxxx_mdout.xtc -maxh 72

         else
		echo "Something went wrong. Exiting the program."
		exit
fi
echo "#####################################################################END################################################################################"
