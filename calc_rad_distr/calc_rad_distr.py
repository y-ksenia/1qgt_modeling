import sys
import string
import os
import numpy as np
from math import pi
import scipy.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm


print('Write the length of DNA in bp:')
bp = int(input().strip())
# print('Write the concentration of NaCl in system:')
# mM = int(input().strip())
mMs = ['0.1','2','5','12']
# print('Write the box size:')
# x_max = int(input().strip())+1
boxes = [1000,1000,1250,950]

# for run in (range(1,2)):
step = 1


# In[34]:
elems = ['DNA','Na','Cl']

for mM, x_max in zip(mMs, boxes):
	print(f'{mM}mM') 
	fname = f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/dsDNA_{bp}bp_{mM}mM_pressure_push.dat'
	frames = np.array([int(i.split('\t')[0].strip()) for i in open(fname)]) // 10000
	for elem in (elems):
	    # print(f'Run_{run}')
	    print(elem)
	    if not os.path.exists(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/rad_distr'):
	        os.mkdir(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/rad_distr')
	    if not os.path.exists(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/rad_distr/{elem}'):
	        os.mkdir(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/rad_distr/{elem}')
	    for frame in tqdm(frames):
	        fname = f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/radius_vectors/{elem}.{frame}.dat'
	        data = np.array([float(i.strip()) for i in open(fname)])

	        x = np.arange(1,x_max+1)
	        pdf = np.array(scipy.stats.gaussian_kde(data).pdf(x))
	        denom = x*x*4*pi
	        pltdt = np.divide(pdf, denom)
	        pltdt /= sum(pltdt)
	        #Сохраняем то, что посчитали
	        file = open(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_from_30nm/rad_distr/{elem}/{frame}.dat', "w")
	        for i,j in enumerate(x):
	            file.write(f'{j}\t{pltdt[i]}\n')
	        file.close()
	        data = None
	        # print(data)

