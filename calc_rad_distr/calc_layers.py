import sys
import string
import os
import numpy as np
from math import pi
import scipy.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm


# print('Write the length of DNA in bp:')
# bp = int(input().strip())
bp=4000
# print('Write the concentration of NaCl in system:')
# mM = float(input().strip())
mM=50
s = 1

# fname = open("{name}/dsDNA_4000bp_mM{mM}_pressure_push.dat") 
# with open(fname, 'r') as f:
#     last = float(f.readlines()[-1]
#              .split('\t')[1])
for run in range(1,6):

    fname = f'../DNA/dsDNA/4000bp/50mM/push_100k_steps/run_{run}/dsDNA_4000bp_50mM_pressure_push.dat'

    rads = np.array([float(i.split('\t')[1].strip()) for i in open(fname)])
    inds =  np.where(rads<=30)[0]
    frames = np.array([int(i.split('\t')[0].strip()) for i in open(fname)])[inds] // 10000
    print(f'Run_{run}\t from {frames[0]} to {frames[-1]} frames')

    num_lay = []
    interlay_dist = []
    #if not os.path.exists(f'../DNA/dsDNA/{bp}bp/{mM}mM/DNA_deriv'):
    #    os.mkdir(f'../DNA/dsDNA/{bp}bp/{mM}mM/DNA_deriv')

    file = open(f'../DNA/dsDNA/{bp}bp/{mM}mM/push_100k_steps/run_{run}_layers.dat', "w")
    # file.write('Frame\t Number of layers\t Average distance between layers\n')

    for frame in tqdm(frames):

        fname = f'../DNA/dsDNA/{bp}bp/{mM}mM/push_100k_steps/rad_distr/run_{run}/DNA/{frame}.dat'
        dens = np.array([float(i.split('\t')[1].strip()) for i in open(fname)])
        x = np.array([float(i.split('\t')[0].strip()) for i in open(fname)]) / 10
        deriv6 = np.zeros(len(dens))
        deriv8 = np.zeros(len(dens))
        
        # file = open(f'dsDNA_4000bp_{mM}mM/DNA_deriv/{frame}.dat', "w")
        for j in range(len(dens)):
            if j > 2 and j < len(dens)-3:
                deriv6[j] = - 1/60 * dens[j-3] + 3/20 * dens[j-2] - 3/4 * dens[j-1] + \
                + 3/4 * dens[j+1] - 3/20 * dens[j+2] + 1/60 * dens[j+3]
            if j > 3 and j < len(dens)-4:
                deriv8[j] = 1/280 * dens[j-4] - 4/105 * dens[j-3] + 1/5 * dens[j-2] - 4/5 * dens[j-1] + \
                + 4/5 * dens[j+1] - 1/5 * dens[j+2] + 4/105 * dens[j+3] - 1/280 * dens[j+4]
            # file.write(str(j) + '\t' + str(deriv6[j]) + '\t' + str(deriv8[j]) + '\n')
        # file.close()
        
        start = int(np.argwhere(deriv6>0)[0])
        finish = int(np.argwhere(dens>0.1e-8)[-1])
        av = np.mean(dens[start:finish])
        count = 0
        x_lay = []

        for i in range(len(dens)):
            if (deriv8[i] > 0 and deriv8[i+s] < 0):# and dens[i] > av:
                count += 1
                x_lay.append(x[i])
        if len(x_lay) > 1:
            dist = []
            for i in range(len(x_lay)-1):
                dist.append(abs(x_lay[i+1]-x_lay[i]))
            av_dist = dist[-1]
            # av.dist = np.mean(dist)
        else:
            av_dist = 0
        av_dist = round(av_dist,2)
        file.write(f'{frame}\t{count}\t{av_dist}\n')
    file.close()
