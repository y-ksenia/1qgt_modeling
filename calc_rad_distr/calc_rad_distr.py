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


step = 1
fname = f'../DNA/dsDNA/{bp}bp/50mM/push_from_30nm/dsDNA_{bp}bp_mM50.0_pressure_push_from_30nm.dat'
with open(fname, 'r') as fh:
    max_frame = int(fh.readlines()[-1].split()[0]) // 10000

# In[34]:

for elem in tqdm(['Na', 'Cl']):
    print(elem)
    if not os.path.exists(f'../DNA/dsDNA/{bp}bp/50mM/{elem}_rad_distr'):
        os.mkdir(f'../DNA/dsDNA/{bp}bp/50mM/{elem}_rad_distr')
    for frame in tqdm(range(0, max_frame, step)):

        fname = f'../DNA/dsDNA/{bp}bp/50mM/radius_vectors/{elem}.{frame}.dat'
        data = np.array([float(i.strip()) for i in open(fname)])

        x = np.arange(10,600)
        pdf = np.array(scipy.stats.gaussian_kde(data).pdf(x))
        denom = x*x*4*pi
        pltdt = np.divide(pdf, denom)

        #Сохраняем то, что посчитали
        file = open(f'../DNA/dsDNA/{bp}bp/50mM/{elem}_rad_distr/{frame}.dat', "w")
        for i,j in enumerate(x):
            file.write(f'{j}\t{pltdt[i]}\n')
        file.close()
        data = None
        # print(data)

