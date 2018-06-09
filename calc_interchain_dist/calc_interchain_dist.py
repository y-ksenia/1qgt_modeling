
# coding: utf-8

# In[25]:

from numpy import *
import math as m

N = 3182
file_name = raw_input()
file = open(file_name + ".xyz", "r")
coords = []
sum = 0
for i, line in enumerate(file):
    parts = line.split()
    if i >= 2:
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        coords.append((x,y,z))
file.close()

#print coords[i][0]
for i in range(N):
    rx = abs(coords[i][0] - coords[i + N][0])
    ry = abs(coords[i][1] - coords[i + N][1])
    rz = abs(coords[i][2] - coords[i + N][2])
    #print rx, ry, rz
    r2 = rx * rx + ry * ry + rz * rz
    r = m.sqrt(r2)
    sum = sum + r
av_r = sum / N
print av_r


# In[ ]:




# In[ ]:



