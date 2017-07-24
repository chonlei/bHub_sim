cd test1
echo "Running test1_1.py: 2 coulped cells test"
python test1_1.py
echo "DONE!"
#echo "Running test1_2.py: coupled cells test"
#python test1_2.py
#echo "DONE!"
cd ../test2
echo "Runnning test2.py: MPI4PY test"
mpiexec -n 16 python test2.py
echo "DONE!"
