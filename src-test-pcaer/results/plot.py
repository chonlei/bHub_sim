import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

files = glob.glob('*.txt')
files.sort()
colour = {
	'1': 'b',
	'2': 'r',
	'3': 'k',
	'4': 'g',
}
for f in files:
    a = np.loadtxt(f)
    a[:, 1] = a[:, 1] / a[0, 1] * 100
    f_id = f[-5]
    if f[-6] == '0':
        plt.plot(a[:, 0], a[:, 1], c=colour[f_id])
    elif f[-6] == '1':
        plt.plot(a[:, 0], a[:, 1], ls='--', c=colour[f_id])
    else:
        print('Not sure what you are doing!')

plt.xlabel('\% inhibition')
plt.ylabel('\% Ca activity')
plt.savefig('all.png')
