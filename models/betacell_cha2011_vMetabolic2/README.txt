To run with python (on Linux)

$ export PYTHONPATH=/usr/loca/nrn/lib/python
if needed.



$ python main.py



OR



$ python

>>> execfile('main.py')

=============================================
SOME NOTES ON THE MODEL:

- main-test.py simulated only two beta cells to check their performances.


- To run with nrngui.exe on Windows, compile the .mod files with the Windows compiler. And use the .mod and .hoc files to build a master.hoc. Also remember to add nrn_load_dll("nrnmech.dll") to the betacell.hoc.


- GLUCOSE LEVEL: to change the glucose level, it should be done in the betacell.hoc file, parameter Glucose. Currently it has been set with glucose level at 8 mM.


- All notations follow the paper Cha et al. 2011 as much as possible.
- Expected result of V[mV]-t[ms] plot of the model is shown in figure_1.png as reference.
- Expected result of Ca[mM]-t[ms] plot of the model is shown in figure_2.png as reference.
- All file naming convention follow Linford's previous model from Jo J 2005. Except morphology.hoc -> betacell.hoc.

