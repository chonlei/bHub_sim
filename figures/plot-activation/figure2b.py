import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

CASE = [1,2,3,4,5]
colour = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']#['r','c','b','k','g']

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
ax.set_xlabel(r"%cells as hub")
ax.set_ylabel(r"[Ca]$_i$ activity [$\mu$M]")

for (idx,subcase) in enumerate(CASE):
    loc = "case%d/"%(subcase)  # location of files
    FILES = glob.glob("%s*.txt"%loc)  # list of files to plot
    for sim in FILES:
        trace = np.loadtxt(sim)
        trace[:,1] = trace[:,1]*1000.
        trace[:,0] += 1./750.
        trace = np.insert(trace, 0, 0, axis=0)
        #trace[:,1] = trace[:,1]/trace[0,1]*100.0  # normalise
        ax.plot(trace[:,0],trace[:,1],c=colour[idx])
    ax.plot(trace[:,0],trace[:,1],c=colour[idx],label='case%d'%subcase)

# add legend
legend = plt.legend()  #plt.legend(["case1", "case2", "case4"])#, loc=4);
frame = legend.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')

plt.savefig("../figures/activation.png",bbox_inches='tight')
plt.savefig("../figures/activation.pdf",format='pdf',bbox_inches='tight')

plt.show()