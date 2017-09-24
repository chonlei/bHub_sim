#!/bin/bash

# $1 : sim folder idx
# $2 : number of subsim folders
# $3 : first Ca .dat file
# $4 : isImitateExp
# $5 : number of simulation batches

pathtofiles="/home/scratch/output_bHub_arcus"
out="outplot"
#mkdir $pathtofiles/sim$1/$out
for subidx in $( seq 0 $2 )
do
	echo plotting subidx $subidx
	# python plotTraces.py $pathtofiles/sim$1/subsim$subidx/$3 $4 $5
	python getS2.py $pathtofiles/sim$1/subsim$subidx/$3 $4 $5 $subidx $1
	#mv $pathtofiles/sim$1/subsim$subidx/image* $pathtofiles/sim$1/$out/sim${subidx}i.png
	#mv $pathtofiles/sim$1/subsim$subidx/whole* $pathtofiles/sim$1/$out/sim${subidx}w.png
done
echo DONE

python quickplot2.py sim$1.txt $2 $pathtofiles/sim$1/subsim$subidx/$3 $1
# rm -rf temp.txt
