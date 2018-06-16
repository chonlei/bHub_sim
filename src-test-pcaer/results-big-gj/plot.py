import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

files = glob.glob('*.txt')
files.sort()
colours = {
	'1': 'b',
	'2': 'r',
	'3': 'k',
	'4': 'g',
}
labels = {
	'1': '-20%',
	'2': '-10%',
	'3': 'CTRL',
	'4': '+10%',
}

for f in files:
    a = np.loadtxt(f)
    a[:, 1] = a[:, 1] / a[0, 1] * 100
    f_id = f[-5]
    if f[-6] == '2':
        plt.plot(a[:, 0], a[:, 1], c=colours[f_id], label=labels[f_id])
    elif f[-6] == '3':
        plt.plot(a[:, 0], a[:, 1], ls='--', c=colours[f_id])
    else:
        print('Not sure what you are doing!')

plt.legend()
plt.xlabel('\% inhibition')
plt.ylabel('\% Ca activity')
plt.savefig('all.png')
