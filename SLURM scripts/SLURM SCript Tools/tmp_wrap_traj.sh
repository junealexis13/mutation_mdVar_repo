#!/bin/bash
#SBATCH --partition=batch
#SBATCH --qos=240c-1h_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name='WrapTrajyyyy'
#SBATCH --output=WrapTraj.%J.out
#SBATCH --error=WrapTraj.%J.err
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

#invoke Rewrap and recenter
if [ -f $WD/xxxx_mdout_center.xtc ]; then
	echo "1 0" | gmx trjconv -s $WD/xxxx_mdout.tpr -f $WD/xxxx_mdout.xtc -o $WD/xxxx_mdout_center.xtc -center -pbc nojump -ur compact
fi

#create an index for ligand RMSD
CMD1="13 & ! a H*"
CMD2="name 21 UNK_Heavy"
CMD3="q"

#bunch of content redirected for input automation =) 
gmx make_ndx -f $WD/xxxx_min.gro -n $WD/xxxx_min.ndx <<< "$CMD1"$'\n'"$CMD2"$'\n'"$CMD3"
gmx rms -s $WD/xxxx_minimization.tpr -f $WD/xxxx_mdout_center.xtc -n $WD/index.ndx -tu ns -o $WD/xxxx_ligandRMSD.xvg <<< "4"$'\n'"21"
