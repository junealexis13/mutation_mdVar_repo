echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# MAIN
WD=$SLURM_SUBMIT_DIR
cd $WD
# PREPARATION PART WAS REDACTED 

if [[ -f $WD/xxxx_mdout.log ]]; then
	gmx mdrun -deffnm xxxx_mdout -nt 6 -ntmpi 2 -ntomp 3 -npme 1 -nb gpu -pme cpu -bonded gpu -s $WD/xxxx_mdout.tpr -x $WD/xxxx_mdout.xtc -cpi xxxx_mdout.cpt -maxh 72
	exit
else
	if [[ -f $WD/nvt_omic.mdp ]] && [[ -f $WD/xxxx_min.gro ]]; then
        	#NVT equilibration
		gmx make_ndx -f $WD/xxxx_min.gro -o $WD/xxxx_min.ndx <<< "1 | 13"$'\n'"16 | 14 | 15"$'\n'"q"
        	gmx grompp -f $WD/nvt_omic.mdp -c $WD/xxxx_min.gro -r $WD/xxxx_min.gro -p $WD/topol.top -n $WD/xxxx_min.ndx -o $WD/xxxx_equi.tpr
        	gmx mdrun -deffnm xxxx_nvt -s $WD/xxxx_equi.tpr -x $WD/xxxx_nvt.xtc
	fi

	if [[ -f $WD/npt_omic.mdp ]] && [[ -f $WD/xxxx_nvt.gro ]]; then
        	#NPT Equilibration
        	gmx grompp -f $WD/npt_omic.mdp -c $WD/xxxx_nvt.gro -r $WD/xxxx_nvt.gro -t $WD/xxxx_nvt.cpt -n $WD/xxxx_min.ndx -o $WD/xxxx_npt.tpr -p $WD/topol.top -maxwarn 3
        	gmx mdrun -deffnm xxxx_npt -s $WD/xxxx_npt.tpr -x $WD/xxxx_npt.xtc

        	#Getting NPT            energy log
        	#echo "18 0" | gmx energy -f $WD/complex_lig_prot_npt.edr -o $WD/complex_lig_prot_pressure.xvg

        	#Getting NPT density log
        	#echo "24 0" | gmx energy -f $WD/complex_lig_prot_npt.edr -o $WD/density.xvg

	fi

	if [[ -f $WD/mdout.mdp ]] && [[ -f xxxx_npt.gro ]]; then
        	#Executing MD RUN
		gmx grompp -f $WD/mdout.mdp -c $WD/xxxx_npt.gro -r $WD/xxxx_npt.gro -t $WD/xxxx_npt.cpt -p $WD/topol.top -n $WD/xxxx_min.ndx -o $WD/xxxx_mdout.tpr
		gmx mdrun -deffnm xxxx_mdout -nt 6 -ntmpi 2 -ntomp 3 -npme 1 -nb gpu -pme cpu -bonded gpu -s $WD/xxxx_mdout.tpr -x $WD/xxxx_mdout.xtc -maxh 72
	fi
fi


echo "#####################################################################END################################################################################"
