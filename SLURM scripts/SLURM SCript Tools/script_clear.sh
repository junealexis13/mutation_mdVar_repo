#!/bin/bash


echo "
SUBMISSION CONTROLLER - Clearing Script
---------------------
Author: June Alexis A. Santos, MSc Student, AdU
NOTE: This script should only be run on front-end.
"


for DIR in $(ls); do
	if [[ -d $(pwd)/$DIR ]]; then
		if ! [[ $DIR == "com_Pr1_C478K" ]]; then 
			cd $(pwd)/$DIR
			rm *.sh
			rm *.mdp
			rm JobName*
			rm *solv*
			rm *box*
			if [[ -f "#topol.top.1#" ]]; then
				mv "#topol.top.1#" "topol.top"
				rm -r #topol.top*
			fi

			cd ..
		fi
	fi
done

echo "DONE"
