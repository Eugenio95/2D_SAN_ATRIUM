#!/bin/bash

## declare an array of strings
declare -a arr=("_HUMAN_noHet5SEP3_Seed1-50s_GRAD_M0_SEPc_F0_14-Nov-2023.txt"
		"_HUMAN_noHet5SEP7_Seed1-50s_GRAD_M0_SEPc_F0_14-Nov-2023.txt"
		"_HUMAN_noHet5SEP9_Seed1-50s_GRAD_M0_SEPc_F0_14-Nov-2023.txt"

	)

count=0

## now loop throughthe above array
for i in "${arr[@]}"
do
	count=$((count+1))
	./2D_parallel_SAN_atrium Sim_param_folder $i &> log$count

done	

