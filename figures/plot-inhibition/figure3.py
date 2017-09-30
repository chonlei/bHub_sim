import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

try:
    CASE = int(sys.argv[1])
except Exception:
    CASE = None
    #raise Exception("Enter case idx to plot as arg")
SUBCASE = ['hub','follower']
colour = ['r','k']

## figure setup
params = {
   'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 10,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
   'figure.figsize': [4.5, 4.5]
   }
plt.rcParams.update(params)

fig = plt.figure(figsize=(6, 4.5))
ax = fig.add_subplot(1, 1, 1)

## axis setup
ax.set_xlabel(r"%hub inhibited")
ax.set_ylabel(r"%activity")


for (idx,subcase) in enumerate(SUBCASE):
    loc = "case%d/%s/"%(CASE,subcase)  # location of files
    FILES = glob.glob("%s*.txt"%loc)  # list of files to plot
    for sim in FILES:
        trace = np.loadtxt(sim)
        trace[:,1] = trace[:,1]/trace[0,1]*100.0  # normalise
        ax.plot(trace[:,0],trace[:,1],c=colour[idx])
    ax.plot(trace[:,0],trace[:,1],c=colour[idx],label='%s'%subcase)

# add legend
legend = plt.legend()  #plt.legend(["case1", "case2", "case4"])#, loc=4);
frame = legend.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')

plt.savefig("../figures/inhibition-%d.png"%CASE,bbox_inches='tight')
plt.savefig("../figures/inhibition-%d.pdf"%CASE,format='pdf',bbox_inches='tight')

plt.show()