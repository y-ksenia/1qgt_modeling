
# coding: utf-8

# In[17]:


import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
from scipy.spatial.distance import cosine
from scipy.linalg import solve
import time
import mdtraj
from tqdm import tqdm
import time
import os


# In[18]:


#Calculate coordinates of a middle point for a base pair as array 
def middle_point(p1, p2):
    return np.mean([p1, p2], 0)

def unit_vector(p1, p2):
    return (p2 - p1) / la.norm((p2 - p1), None)

def calculate_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    angle = np.arccos(np.dot(v1,v2)/(la.norm(v1, None)*LA.norm(v2, None)))
    return angle

def сalculate_central_line(dna):
    central_points = []
    N = len(dna)//2
    
    for i in range(N):
        #The matrix of normals
        A = np.zeros((3,3))
        
        #The vector - right side of system
        b = np.zeros(3)
        
        #Get the positions of the i-pair
        i1 = i
        i2 = i + N
        p1 = dna.positions[i1]    
        p2 = dna.positions[i2]
        
        #Calculate the middle point of the pair and the normalized vector between it's point
        k1 = middle_point(p1, p2)
        n1 = unit_vector(p1,p2)
        A[0] = n1
        b[0] = np.dot(n1,k1)
        
        #Find the next pair - it has bigger angle with the i-pair
        max_cos = 0.
        j_max = 0
        n = 0
        for j in range(-3, 3):
            if (j == 0) or (i + j < 0) or (i + j >= N):
                continue
                
            i3 = i + j
            i4 = i + j + N
            p3 = dna.positions[i3]    
            p4 = dna.positions[i4]
            
            n2 = unit_vector(p3,p4)
            
            if cosine(n1, n2) > max_cos:
                max_cos = cosine(n1, n2)
                j_max = j
            
            if cosine(n1, n2) > 0.05:
                k2 = middle_point(p3, p4)
                A[1] = n2
                b[1] = np.dot(n2,k2)
                
                n3 = np.cross(n1,n2)
                A[2] = n3
                b[2] = np.dot(n3,k1)

                c_point = solve(A,b)
                n += 1
                
        if n == 0:
            i3 = i + j_max
            i4 = i + j_max + N
            p3 = dna.positions[i3]    
            p4 = dna.positions[i4]
            
            n2 = unit_vector(p3,p4)
            k2 = middle_point(p3, p4)
            A[1] = n2
            b[1] = np.dot(n2,k2)
            
            n3 = np.cross(n1,n2)
            A[2] = n3
            b[2] = np.dot(n3,k1)
            
            c_point = solve(A,b)
            n += 1
            
        central_points.append(c_point)
    return central_points


# In[26]:

#print('Write the length of DNA in bp:')
bp = 4000
#print('Write the concentration of NaCl in system:')
#conc = int(input().strip())
cc = [0.1,2,5,12]

for conc in cc:	

	print(f'System: {bp} bp, {conc} mM')
#	print(f'Run {run} \n')
	path = f'../DNA/dsDNA/{bp}bp/{conc}mM/push_from_30nm'
	if not os.path.exists(f'{path}/central_line'):
		os.mkdir(f'{path}/central_line')
	nameIN = f'dsDNA_{bp}bp_{conc}mM_push'
	nameOUT = f'central_line/dsDNA_{bp}bp_{conc}mM'
	u = MDAnalysis.Universe(f'{path}/{nameIN}.psf', 
	                        f'{path}/{nameIN}.dcd')
	print('Done reading psf and dcd...')
	dna = u.select_atoms('resname DNA')


	# In[27]:


	#Write xyz-file for central line
	print('Start writing xyz...')
	for_xyz = np.array(сalculate_central_line(dna), dtype = np.float32)
	xyz_output = open(f'{path}/{nameOUT}.xyz', "w")
	xyz_output.write(str(len(for_xyz)) + "\n")
	xyz_output.write("Generate DNA central line"+"\n")
	for ia in for_xyz:
	    xyz_output.write("C\t")	
	    xyz_output.write(str(round(ia[0],7))+"\t")
	    xyz_output.write(str(round(ia[1],7))+"\t")
	    xyz_output.write(str(round(ia[2],7))+"\n")
	xyz_output.close()
	print('Done writing xyz...')
	#Запись углов в dat файл
	#angles_output = open(sys.argv[3], "w")
	#for ang in angles:
	#	angles_output.write(str(ang)+"\n")
	#angles_output.close()


	# In[ ]:

	print('Start writing dcd...')
	f = mdtraj.formats.DCDTrajectoryFile(f'{path}/{nameOUT}.dcd', 'w')
	for ts in tqdm(u.trajectory):
	    # start_time = time.time()

	    #Comment next two lines if there is no need in PBC
	    # dna.positions = dna.positions - dna.positions[50]
	    # dna.positions = dna.positions - ((dna.positions + box_lenght / 2) // box_lenght) * box_lenght

	    c_coords = np.array(сalculate_central_line(dna), dtype = np.float32)
	    f.write(c_coords)
	    # time.sleep(0.01)
	f.close()
	print('Done writing dcd...')

