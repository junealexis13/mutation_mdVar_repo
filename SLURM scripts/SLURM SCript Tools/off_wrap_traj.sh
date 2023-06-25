#!/bin/bash
#create an index for ligand RMSD
CMD1="13 & ! a H*"
CMD2="name 21 UNK_Heavy"
CMD3="q"
CMD4="r 436-506"

module load gromacs/2021

if [ -f "$WD/alpha_WT_complex_mdout.log" ]; then
	if grep -Fq "Writing checkpoint, step 50000000" $WD/alpha_WT_complex_mdout.log; then
		#invoke Rewrap and recenter
		if [ ! -f $WD/alpha_WT_complex_mdout_center.xtc ]; then
			echo "1 0" | gmx trjconv -s $WD/alpha_WT_complex_mdout.tpr -f $WD/alpha_WT_complex_mdout.xtc -o $WD/alpha_WT_complex_mdout_center.xtc -center -pbc nojump -ur compact
		fi


		if [ -f $WD/alpha_WT_complex_mdout_center.xtc ]; then
			#compute RMSD using bunch of content redirected for input automation =) 
			echo "Found CENTERED FILES. PROCEEDING WITH DATA EXTRACTION

			gmx make_ndx -f $WD/alpha_WT_complex_min.gro -n $WD/alpha_WT_complex_min.ndx <<< "$CMD1"$'\n'"$CMD2"$'\n'"$CMD3"
			gmx rms -s $WD/alpha_WT_complex_minimization.tpr -f $WD/alpha_WT_complex_mdout_center.xtc -n $WD/alpha_WT_complex_min.ndx -tu ns -o $WD/alpha_WT_complex_ligandRMSD.xvg <<< "4"$'\n'"21"

			#computing the RMSF
			gmx make_ndx -f $WD/alpha_WT_complex_min.gro -n $WD/alpha_WT_complex_min.ndx <<< "$CMD4"$'\n'"$CMD3"
			#C-ALPHA RMSF
			gmx rmsf -s $WD/alpha_WT_complex_minimization.tpr -f $WD/alpha_WT_complex_mdout_center.xtc -n $WD/alpha_WT_complex_min.ndx -tu ns -o $WD/alpha_WT_complex_ligandRMSF.xvg <<< "3"
			#RBM RMSF
			gmx rmsf -s $WD/alpha_WT_complex_minimization.tpr -f $WD/alpha_WT_complex_mdout_center.xtc -n $WD/alpha_WT_complex_min.ndx -tu ns -o $WD/alpha_WT_complex_ligandRMSF.xvg <<< "22"
			#GET SASA
			gmx sasa -f alpha_WT_complex_mdout_center.xtc -s alpha_WT_complex_minimization.tpr -n $WD/alpha_WT_complex_min.ndx -o alpha_WT_complex_sasa.xvg -odg alpha_WT_complex_solvation_free_energy.xvg -or alpha_WT_complex_avg_perres_sasa.xvg -tv alpha_WT_complex_total_v.xvg <<< "19"

		fi
	else
		echo "Simulation not yet finished"
		exit
	fi
fi