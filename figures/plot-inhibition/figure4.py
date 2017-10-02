import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

try:
    CASE = int(sys.argv[1])
except Exception:
    CASE = 5
    #raise Exception("Enter case idx to plot as arg")
SUBCASE = ['hub','follower']
colour = ['tab:red','tab:grey']

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
ax.set_xlim([-5, 125])
ax.set_ylim([-5, 105])


## function to get quartiles
def perc(data):
    data = data.T  # coomment this out if needed
    median = np.zeros(data.shape[1])
    perc_25 = np.zeros(data.shape[1])
    perc_75 = np.zeros(data.shape[1])
    for i in range(0, len(median)):
        median[i] = np.median(data[:, i])
        perc_25[i] = np.percentile(data[:, i], 25)
        perc_75[i] = np.percentile(data[:, i], 75)
    return median, perc_25, perc_75


for (idx,subcase) in enumerate(SUBCASE):
    loc = "case%d/%s/"%(CASE,subcase)  # location of files
    FILES = glob.glob("%s*.txt"%loc)  # list of files to plot
    alltraces = []
    for sim in FILES:
        trace = np.loadtxt(sim)
        trace[:,1] = trace[:,1]/trace[0,1]*100.0  # normalise
        #ax.plot(trace[:,0],trace[:,1],c=colour[idx])
        alltraces.append(trace[:31,1])
    alltraces = np.array(alltraces)
    meantr = np.mean(alltraces,0)
    stdtr = np.std(alltraces,0)
    ax.plot(trace[:31,0],meantr,c=colour[idx],label='%s'%subcase)
    #med_Inorm, perc_25_Inorm, perc_75_Inorm = perc(alltraces.T)
    ax.fill_between(trace[:31,0], meantr-stdtr, meantr+stdtr, alpha=0.25, lw=0, color=colour[idx])
    #ax.fill_between(trace[:31,0], perc_25_Inorm, perc_75_Inorm, alpha=0.25, lw=0, color=colour[idx])

# add legend
legend = plt.legend()  #plt.legend(["case1", "case2", "case4"])#, loc=4);
frame = legend.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')

plt.savefig("../figures/inhibition-%d_av.png"%CASE,bbox_inches='tight')
plt.savefig("../figures/inhibition-%d_av.pdf"%CASE,format='pdf',bbox_inches='tight')

plt.show()