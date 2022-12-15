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
"


for DIR in $(ls); do
	if [[ $DIR != "com_Pr1_C478K" ]]; then 
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
			else
				echo "Script was not submitted. Maybe this is because of inappropriate flag or the files doesn't exist."
			fi
			cd .. 
		fi
	fi
done



