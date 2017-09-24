#!/bin/bash

# create output directory
outputdir=$(python outputSetup.py ../output/ >&1)

# run main simulations
nsim=7
for i in {0..4}; do
	for (( j=0; j<$nsim; j++ )); do
		nohup python main.py ${outputdir} 3 2 1 $(( $nsim*$i + $j )) &
		sleep 2s
	done
	wait $!
done

