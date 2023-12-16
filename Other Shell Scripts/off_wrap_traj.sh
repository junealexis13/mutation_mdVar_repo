#!/bin/bash
#create an index for ligand RMSD
CMD1="13 & ! a H*"
CMD2="name 21 UNK_Heavy"
CMD3="q"
CMD4="r 436-506"

if [ -f "./com_Pr1_E484A_mdout.log" ]; then
	if grep -Fq "Writing checkpoint, step 50000000" ./com_Pr1_E484A_mdout.log; then
		#invoke Rewrap and recenter
		if [ ! -f ./com_Pr1_E484A_mdout_center.xtc ]; then
			echo "1 0" | gmx trjconv -s ./com_Pr1_E484A_mdout.tpr -f ./com_Pr1_E484A_mdout.xtc -o ./com_Pr1_E484A_mdout_center.xtc -center -pbc nojump -ur compact
		fi


		if [ -f ./com_Pr1_E484A_mdout_center.xtc ]; then
			#compute RMSD using bunch of content redirected for input automation =) 
			echo "Found CENTERED FILES. PROCEEDING WITH DATA EXTRACTION"

			gmx make_ndx -f ./com_Pr1_E484A_min.gro -n ./com_Pr1_E484A_min.ndx <<< "$CMD1"$'\n'"$CMD2"$'\n'"$CMD3"
			gmx rms -s ./com_Pr1_E484A_minimization.tpr -f ./com_Pr1_E484A_mdout_center.xtc -n ./index.ndx -o ./com_Pr1_E484A_ligandRMSD.xvg <<< "4"$'\n'"21"

			#computing the RMSF
			gmx make_ndx -f ./com_Pr1_E484A_min.gro -n ./index.ndx <<< "$CMD4"$'\n'"$CMD3"
			#C-ALPHA RMSF
			gmx rmsf -s ./com_Pr1_E484A_minimization.tpr -f ./com_Pr1_E484A_mdout_center.xtc -n ./index.ndx  -o ./com_Pr1_E484A_cAlphaRMSF.xvg -res <<< "3"
			#RBM RMSF
			gmx rmsf -s ./com_Pr1_E484A_minimization.tpr -f ./com_Pr1_E484A_mdout_center.xtc -n ./index.ndx -o ./com_Pr1_E484A_ligandRMSF.xvg -res <<< "22"
			#GET SASA
			gmx sasa -f com_Pr1_E484A_mdout_center.xtc -s com_Pr1_E484A_minimization.tpr -n ./index.ndx -o com_Pr1_E484A_sasa.xvg -odg com_Pr1_E484A_solvation_free_energy.xvg -or com_Pr1_E484A_avg_perres_sasa.xvg -tv com_Pr1_E484A_total_v.xvg <<< "19"

			#Calculating Complex Gyration
			if ! [ -f *ProtLig_gyrate* ]; then
				echo "19" | gmx gyrate -f com_Pr1_E484A_mdout_center.xtc -s com_Pr1_E484A_minimization.tpr -n index.ndx -o com_Pr1_E484A_ProtLig_gyrate.xvg
			else
				echo "Gyration Data has been found. Skipping this step..."
			fi

			#Calculating Protein Gyration
			if ! [ -f *ProtOnly_gyrate* ]; then
				echo "1" | gmx gyrate -f com_Pr1_E484A_mdout_center.xtc -s com_Pr1_E484A_minimization.tpr -n index.ndx -o com_Pr1_E484A_ProtOnly_gyrate.xvg
			else
				echo "Gyration Data has been found. Skipping this step..."
			fi

			#get the hydrogen bond distribution
			if [[ ! -f hbnum.xvg ]] && [[ ! -f hbdist.xvg ]]; then
				echo "1 13" | gmx hbond -f com_Pr1_E484A_mdout_center.xtc -s com_Pr1_E484A_minimization.tpr -n index.ndx -num -dist
				mv hbnum.xvg com_Pr1_E484A_hbnum.xvg
				mv hbdist.xvg com_Pr1_E484A_hbdist.xvg
			else
				echo "Hydrogen Bond Distribution data has been found. Skipping this step..."
			fi


		fi
	else
		echo "Simulation not yet finished"
		exit
	fi
fi
