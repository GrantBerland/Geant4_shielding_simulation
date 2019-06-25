#!/bin/bash

read -p "Enter number of simulations to run: " numberSims
read -p "Number of cores to use: " cores
echo

loops=$((numberSims / cores));

echo "Running $numberSims simulations in $loops loops..."

for sims in $(seq 1 $loops); 

	do parallel ./parallel_runfile ::: {1..$cores} >> output.log; 
	echo "Simulation $((sims*cores))/$numberSims complete"; 

done;

echo "All simulations complete"


