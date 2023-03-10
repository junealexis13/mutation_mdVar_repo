#!/bin/bash
echo "
SUBMISSION CONTROLLER
---------------------
Author: June Alexis A. Santos, MSc Student, AdU
NOTE: This script should only be run on front-end.
"

gpu-mode () {
echo "#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --qos=12c-1h_2gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --exclude=saliksik-gpu-02
#SBATCH --job-name='yyyy'
#SBATCH --output=JobName.%J.out
#SBATCH --error=JobName.%J.err
#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
#SBATCH --mail-type=ALL
module load gromacs/2021.4_cuda-11.4.3

"
}

cpu-mode () {
echo "#!/bin/bash
#SBATCH --partition=batch
#SBATCH --qos=240c-1h_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name='yyyy'
#SBATCH --output=JobName.%J.out
#SBATCH --error=JobName.%J.err
#SBATCH --mail-user=june.alexis.santos@adamson.edu.ph
#SBATCH --mail-type=ALL
module load gromacs/2021

"
}

WD=$(pwd)

#Preemptively sets the workmode. Change the value based on your liking
declare -A jobHash=( ["Pr1"]="gpu" ["Pr2"]="cpu" ["Pr3"]="cpu" ["Pr4"]="cpu" ["Pr5"]="cpu" ["Pr6"]="cpu" ["Pr7"]="cpu" ["Pr8"]="cpu" )

echo "
Make sure to check the script before SBATCH Submission. =)
"

for JOBS in "${!jobHash[@]}"; do
	for DIRS in $(ls); do
		if ! [[ $DIRS == "com_Pr1_C478K" ]]; then
			if [[ -d ./$DIRS ]] && [[ $DIRS == *"$JOBS"* ]]; then
				if [[ "${jobHash[$JOBS]}" == "cpu" ]]; then
					echo "RUNNING FOR $DIRS SCRIPT CONTROL ALGORITHM - CPU" 
					cd ./$DIRS/ && echo "copying files..." && echo "Doing $DIRS"
					echo "Copying necessary templat files..."
					cp ../*.mdp .
					echo "Creating slurm submit script...for $DIRS"
					cpu-mode > omic_submit_pt1.sh
					cpu-mode > omic_submit_pt2.sh
					echo "Done!"
					chmod +x *sh
					echo "Appending code chunks..."
					cat ../tmp_omic_run_complex.sh >> omic_submit_pt1.sh
					cat ../tmp_omic_run_complex_cpu_pt2.sh >> omic_submit_pt2.sh
					echo "applying in line script modifications..."
					sed -i "s/xxxx/$DIRS/g" omic_submit_pt1.sh
					sed -i "s/yyyy/$JOBS/g" omic_submit_pt1.sh
					sed -i "s/xxxx/$DIRS/g" omic_submit_pt2.sh
					sed -i "s/yyyy/$JOBS/g" omic_submit_pt2.sh 
					cd ..
				elif [[ "${jobHash[$JOBS]}" == "gpu" ]]; then
                                	echo "RUNNING FOR $DIRS SCRIPT CONTROL ALGORITHM - GPU"
					echo "Creating slurm submit script... for $DIRS"
					cd ./$DIRS/ && echo "copying files..."
                                	gpu-mode > omic_submit_pt1.sh
                                	gpu-mode > omic_submit_pt2.sh
					echo "Copying necessary files..."
					cp ../*.mdp .
                                	chmod +x *sh
					echo "Appending code chunks..."
                                	cat ../tmp_omic_run_complex.sh >> omic_submit_pt1.sh
                                	cat ../tmp_omic_run_complex_gpu_pt2.sh >> omic_submit_pt2.sh
					echo "applying in line script modifications..."
					sed -i "s/xxxx/$DIRS/g" omic_submit_pt1.sh
					sed -i "s/yyyy/$JOBS/g" omic_submit_pt1.sh
					sed -i "s/xxxx/$DIRS/g" omic_submit_pt2.sh
					sed -i "s/yyyy/$JOBS/g" omic_submit_pt2.sh 
                                	cd ..
					echo "Done!"
				else
					continue
				fi
			fi
		fi
	done
done



echo "JOB FINISHED!"
