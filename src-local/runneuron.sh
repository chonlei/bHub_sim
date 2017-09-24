#!/bin/bash

# create output directory
outputdir=$(python outputSetup.py ../output/ >&1)

# run main simulations
nsim=6
for i in {0..4..1};
do for j in {0..${nsim}..1}; do python main.py ${outputdir} 3 2 1 $(( ($nsim+1)*$i + $j )); done; wait $!;
done;

