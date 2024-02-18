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


# PREPARATION PART WAS REDACTED 

if [[ -f $WD/xxxx.gro ]]; then
        echo "Found $WD/xxxx.gro"
        #configuration for solvation
        gmx editconf -f $WD/xxxx_complex.gro -o $WD/xxxx_box.gro -d 1.0 -bt cubic
        gmx solvate -cp $WD/xxxx_box.gro -cs spc216.gro -o $WD/xxxx_solv.gro  -p $WD/topol.top
fi

if [[ -f $WD/xxxx_solv.gro ]] && [[ -f $WD/ions.mdp ]] ; then
	#grompping for GENION
	gmx grompp -f $WD/ions.mdp -c $WD/xxxx_solv.gro -o $WD/xxxx_solv_ions.tpr -p $WD/topol.top -maxwarn 3
        echo "15" | gmx genion -s $WD/xxxx_solv_ions.tpr -o $WD/xxxx_solv_ions.gro -p $WD/topol.top -pname SOD -nname CLA -neutral -conc 0.15
fi

if [[ -f $WD/minimization.mdp ]] && [[ -f $WD/xxxx_solv_ions.gro ]]; then
        #energy minimization
        gmx grompp -f $WD/minimization.mdp -c $WD/xxxx_solv_ions.gro -p $WD/topol.top -o $WD/xxxx_minimization.tpr
        gmx mdrun -ntmpi 16 -ntomp 5 -npme 4 -s $WD/xxxx_minimization.tpr -deffnm $WD/xxxx_min -x $WD/xxxx_minimization.xtc

        #getting the minimization energy
        #echo "10 0" | gmx energy -f $WD/xxxx_min.edr -o $WD/xxxx_min.xvg
		
fi


exit 
echo "#####################################################################END################################################################################"
