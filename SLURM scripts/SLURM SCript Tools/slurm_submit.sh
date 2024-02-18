#!/bin/bash
PT=$1

echo "
SLURM SUBMISSION SCRIPT (PT1)
-----------------------------
Author: June Alexis A. Santos, MSc Student, AdU

Note: This script runs on front-end
Required flags:
Scripts to run:
PT1 --> Run only part 1 script to slurm submit
PT2 --> Run only part 2 script to slurm submit
PT3 --> Rerun all CPT Scripts (Checkpoint rerun)
PT4 --> Run Trajectory Wrapping
PT5 --> Run Trajectory Analysis
"


for DIR in $(ls); do
	if [[ -d $(pwd)/$DIR ]] ; then
		cd $(pwd)/$DIR/
		if [[ $PT == "pt1" ]] || [[ $PT == "PT1" ]]; then 
			echo "--- Submitting $DIR/omic_submit_pt1.sh ---" 
			sbatch omic_submit_pt1.sh
			echo "Done!"
			sleep 3
		elif [[ $PT == "pt2" ]] || [[ $PT == "PT2" ]]; then
			echo "--- Submitting $DIR/omic_submit_pt2.sh ---"
			sbatch omic_submit_pt2.sh
			echo "Done!"
			sleep 3
		elif [[ $PT == "pt3" ]] || [[ $PT == "PT3" ]]; then
			echo "--- Submitting $DIR/continue_from_cpy.sh ---"
			sbatch continue_from_cpy.sh 
			echo "Done!"
			sleep 3
		elif [[ $PT == "pt4" ]] || [[ $PT == "PT4" ]]; then
			if [ -f wrap_traj.sh ]; then
				echo "--- Submitting $DIR/wrap_traj.sh ---"
				sbatch wrap_traj.sh 
				echo "Done!"
				sleep 3
			else
				echo "Wrapping script not found. Skipping..."
			fi
		elif [[ $PT == "pt5" ]] || [[ $PT == "PT5" ]]; then
			if [ -f data_extraction.sh ]; then
				echo "--- Submitting $DIR/data_extraction.sh ---"
				sbatch data_extraction.sh 
				echo "Done!"
			else
				echo "Wrapping script not found. Skipping..."
			fi
		else
			echo "Script was not submitted. Maybe this is because of inappropriate flag or the files doesn't exist."
		fi
		cd .. 
		squeue -u junealexis.santos
	fi
done



