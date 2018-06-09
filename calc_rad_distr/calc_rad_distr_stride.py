import sys
import string
import os
import numpy as np
from math import pi
import scipy.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm

# In[30]:


# In[31]:
print('Write the system name')
name = input().strip()

#Указываем: шаг по фреймам, максимальное количество фрейм (как автоматизировать?)
stride=20
step = 1
fname = f'../DNA_packaging/dsDNA/mM50_length/10000bp/push_from_30nm/dsDNA_10000bp_mM50.0_pressure_push_from_30nm.dat'
with open(fname, 'r') as fh:
    max_frame = int(fh.readlines()[-1].split()[0]) // 10000

# In[34]:

data = None
for elem in (['DNA']):
    print(elem)
    if not os.path.exists(f'{name}/{elem}_rad_distr_stride{sride}'):
        os.mkdir(f'{name}/{elem}_rad_distr_stride{sride}')
    for frame in tqdm(range(0, max_frame, step)):
        #Считываем координаты всех частиц
        fname = f'{name}/distr/{elem}.{frame}.dat'
        if data is None:
            data = np.array([float(i.strip()) for i in open(fname)])
        else:
            data += np.array([float(i.strip()) for i in open(fname)])

        x = np.arange(10,600)
        if frame%stride == 0:
            data /= stride
            pdf = np.array(scipy.stats.gaussian_kde(data).pdf(x))
            denom = x*x*4*pi
            pltdt = np.divide(pdf, denom)

            #Сохраняем то, что посчитали
            file = open(f'{name}/{elem}_rad_distr_stride{sride}/{frame}.dat', "w")
            for i,j in enumerate(x):
                file.write(f'{j}\t{pltdt[i]}\n')
            file.close()
            data = None
            # print(data)

