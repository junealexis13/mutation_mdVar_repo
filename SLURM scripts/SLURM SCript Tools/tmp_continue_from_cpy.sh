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



if [ -f "$WD/xxxx_mdout.log" ]; then
	if grep -Fq "Writing checkpoint, step 50000000" $WD/xxxx_mdout.log; then
		exit
	else
		gmx mdrun -deffnm xxxx_mdout -nb gpu -pme cpu -bonded gpu -s $WD/xxxx_mdout.tpr -x $WD/xxxx_mdout.xtc  -cpi $WD/xxxx_mdout.cpt -maxh 72
		sbatch continue_from_cpy.sh
		exit
	fi
fi
	
