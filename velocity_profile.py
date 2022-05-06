#import necessary modules
import sys
import numpy as np
import tecplot_io as tec
import os
import pandas as pd
#change working directory
path = 'd:\post\Project\Fixed_Turbine'
os.chdir(path)
#necessary variables
r = 40.0
D = 2.0*r
nx = 192
ny = 192
nz = 64
nvar = 6
nturbine = 16
dt = 0.68543297937
tis = 200
tie = 15000
tii = 200
nti = int((tie - tis) / tii + 1)
tii = 200
U_hub = 9.45
#frequency of wave
omega = 0.0498561155567
#wave period
T = 2*np.pi/omega
#turbine hub height
H_hub = 70.0
#Rotational angular period
T_turb = 42.84
U_star = 0.356
H_hub = 70
#infinite velocity
U_i = 11.5258407161
#input folder
foldername =  "./POST_UI_2D1_0001/"
#output folder
outputfolder = "./" + "_wake/"
#declare necessary arrays
U_mean = np.zeros((nz))
W_mean = np.zeros((nz))
U_prime_square = np.zeros((nz))
data2 = np.zeros((nz,2))
data3 = np.zeros((nz,2))
data4 = np.zeros((nz,2))
#------------------------------------------------------------------------------#
#clean data
#get mean compponent

for it in range(nti):
	ti = tis + tii * it
	filename = foldername + "POST_UI_2D1_" + "{:010d}".format(ti) + "_0001.DAT"
	data0 = tec.tecplot_reader(filename, [nz, ny, nvar], 2)
	data0 = data0.reshape([nz,ny,nvar])
	X = data0[:,:,0]
	Y = data0[:,:,1]
	Z = data0[:,:,2]
	U = data0[:,:,3]
	V = data0[:,:,4]
	W = data0[:,:,5]
#  uncomment if you want to use mean velocity horizontal average
#	U_mean = np.mean(U, axis=1)
	Z = np.mean(Z, axis=1)
#  get velocity profile at center cross section
	U_cross = U[:,96]
	W_cross = W[:,96]
	U_mean = U_mean + U_cross
	W_mean = W_mean + W_cross
#time average velocity profile
U_mean = U_mean / nti
W_mean = W_mean / nti

#checked
#print(U_mean)
#print (U_mean_all)
#print(U[1,:])
#print(Z)

#get fluctuating component
for it in range(nti):
	ti = tis + tii * it
	filename = foldername + "POST_UI_2D1_" + "{:010d}".format(ti) + "_0001.DAT"
	data1 = tec.tecplot_reader(filename, [nz, ny, nvar], 2)
	data1 = data1.reshape([nz,ny,nvar])
	X = data1[:,:,0]
	Y = data1[:,:,1]
	U = data1[:,:,3]
	V = data1[:,:,4]
	W = data1[:,:,5]
#  get velocity profile at center cross section
	U = U[:,96]
	W = W[:,96]
	U_fluc = (U-U_mean)
	W_fluc = (W-W_mean)
	Reynolds_stress = - U_fluc*W_fluc
	U_prime_square = U_prime_square + U_fluc**2


U_prime_rms = np.sqrt(U_prime_square / nti)
Reynolds_stress_tm = Reynolds_stress / nti


#normalize data
#save data
data2[:,0] = U_mean*(U_i/U_hub)
data2[:,1] = Z/D
data3[:,0] = U_prime_rms*(U_i/U_star)
data3[:,1] = Z/D
data4[:,0] = Reynolds_stress_tm/(U_star**2)
data4[:,1] = Z/D

#save data into tecplot format
f0 = open( outputfolder + "Vertical profile time-average velocity.plt",'w')
f0.write("VARIABLES = u/U_hub, z/D  \n")
f1 = open( outputfolder + "root mean square velocity.plt",'w')
f1.write("VARIABLES = u_rms/U_hub, z/D  \n")
f2 = open( outputfolder + "Reynolds stress.plt",'w')
f2.write("VARIABLES = (-uv), z/D  \n")
np.savetxt(f0, data2)
np.savetxt(f1, data3)
np.savetxt(f2, data4)
f0.close()
f1.close()
f2.close()