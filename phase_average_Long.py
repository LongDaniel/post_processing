import numpy as np
import os
import tecplot_io as tec
import h5py


def diff_central(x, y):
  x0 = x[:-2]
  x1 = x[1:-1]
  x2 = x[2:]
  y0 = y[:-2]
  y1 = y[1:-1]
  y2 = y[2:]
  f = (x2 - x1)/(x2 - x0)
  f1 = (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
  f2 = x.copy()
  f2[1:-1] = f1
  f2[0] = f1[0]
  f2[-1] = f1[-1]
  return f2

def find_index(_z, _limits):
  _n = len(_z)
  _i_min = 0
  _i_max = _n - 1
  _limits2 = np.zeros(2)
  if isinstance(_limits, float):
    _limits2[0:2] = _limits
  else:
    _limits2 = _limits
          
  for i in range(_n):
    if _z[i]<_limits2[0] and i>_i_min :
      _i_min = i
    if _z[i]>_limits2[1] and i<_i_max :
      _i_max = i
  #print('zlimits='+str(_limits))
  #print('i_min='+str(_i_min)+', i_max='+str(_i_max))
  
  if isinstance(_limits, float):
    return _i_min
  else:
    return _i_min, _i_max

#function calculate the phase average
#input: folder name, phi, serial number, start time, end time, time interval
#output: phase average for that folder
def phase_average(folder, phi, t, tis, tie, tii):
    #change working directory
    os.chdir(folder)
    print("current working directory is: {0}".format(os.getcwd()))
    #read data
    x_phase = np.zeros([NPZ, 48])
    z_phase = np.zeros([NPZ, 48])
    u_phase = np.zeros([NPZ, 48])
    v_phase = np.zeros([NPZ, 48])
    w_phase = np.zeros([NPZ, 48])
    uu_phase = np.zeros([NPZ, 48])
    vv_phase = np.zeros([NPZ, 48])
    ww_phase = np.zeros([NPZ, 48])
    uw_phase = np.zeros([NPZ, 48])
    nti = int((tie - tis) / tii + 1)
    for it in range (nti):
        ti = tis + tii * it
        fname = 'Phase_Average_Phi{:01d}_{:04d}_{:08d}.dat'.format(phi,t,ti)
        f = tec.tecplot_reader(fname, [NPZ, NPX, 6], 2)
        f = f.reshape([NPZ, NPX, 6])
        x = f[:,:,0]
        z = f[:,:,2]
        u = f[:,:,3]
        v = f[:,:,4]
        w = f[:,:,5]
        x_phase = x_phase + x[:,48:96]
        z_phase = z_phase +z[:,48:96]
        u_phase = u_phase + u[:,48:96]
        v_phase = v_phase + v[:,48:96]
        w_phase = w_phase + w[:,48:96]
        #uu_phase = uu_phase + u_phase**2
        #vv_phase = vv_phase +  v_phase**2
        #ww_phase =  ww_phase + w_phase**2
        #uw_phase =  uw_phase + u_phase*w_phase
    #time average
    x_phase = x_phase / int((tie - tis) / tii + 1)
    z_phase = z_phase / int((tie - tis) / tii + 1)
    u_phase = u_phase / int((tie - tis) / tii + 1)
    v_phase = v_phase / int((tie - tis) / tii + 1)
    w_phase = w_phase / int((tie - tis) / tii + 1)
    #uu_phase = uu_phase / int((tie - tis) / tii + 1)
    #vv_phase = vv_phase / int((tie - tis) / tii + 1)
    #ww_phase = ww_phase / int((tie - tis) / tii + 1)
    #uw_phase = uw_phase / int((tie - tis) / tii + 1)
    for it in range (nti):
        ti = tis + tii * it
        fname = 'Phase_Average_Phi{:01d}_{:04d}_{:08d}.dat'.format(phi,t,ti)
        f = tec.tecplot_reader(fname, [NPZ, NPX, 6], 2)
        f = f.reshape([NPZ, NPX, 6])
        x = f[:,:,0]
        z = f[:,:,2]
        u = f[:,:,3]
        v = f[:,:,4]
        w = f[:,:,5]
        uu_phase = uu_phase + (u[:,48:96] - u_phase)**2
        vv_phase = vv_phase + (v[:,48:96] - v_phase)**2
        ww_phase = ww_phase + (w[:,48:96] - w_phase)**2
        uw_phase = uw_phase + (u[:,48:96] - u_phase)*(w[:,48:96] - w_phase)

    uu_phase = uu_phase / nti
    vv_phase = vv_phase / nti
    ww_phase = ww_phase / nti
    uw_phase = uw_phase / nti

    #reshape data
    x_phase = x_phase.reshape([NPZ*48])
    z_phase = z_phase.reshape([NPZ*48])
    u_phase = u_phase.reshape([NPZ*48])
    v_phase = v_phase.reshape([NPZ*48])
    w_phase = w_phase.reshape([NPZ*48])
    uu_phase = uu_phase.reshape([NPZ*48]) 
    vv_phase = vv_phase.reshape([NPZ*48]) 
    ww_phase = ww_phase.reshape([NPZ*48]) 
    uw_phase = uw_phase.reshape([NPZ*48]) 
    #store data into array
    data = np.zeros([NPZ*48, 9])
    data[:,0] = x_phase
    data[:,1] = z_phase
    data[:,2] = u_phase*U/U_star
    data[:,3] = v_phase*U/U_star
    data[:,4] = w_phase*U/U_star
    data[:,5] = uu_phase*((U/U_star)**2)
    data[:,6] = vv_phase*((U/U_star)**2)
    data[:,7] = ww_phase*((U/U_star)**2)
    data[:,8] = uw_phase*((U/U_star)**2)
    #save data
    return data

#parameters
#Declare some variable 
nturbinex = 4
nturbiney = 4
nturbine = 16
#Diameter of turbine
D = 80
#Other variable
Sx = 7
Sy = 7
#kappa = 0.4
#nu = 1.511e-5
#PEX = 1.45444104333
#PEY = 8.72664625997
#hbar = 0.46
#uinfty = 2.54390548295
dt = 0.2941315855
#Rotational angular period
T_turb = 42.84
U_star = 0.356
H_hub = 70
#Mean finite velocity
U = 11.5258407161
tis = 100
tie = 15000
tii = 100
nti = int((tie - tis) / tii + 1)
NPX = 192
NPY = 192
NPZ = 65
path = path = 'd:\post\Fixed_new'
os.chdir(path)

#process phase average
data1 = phase_average(path+'\PhaseAveragePhi0_0001', 0, 1, 0, 900, 3)
data2 = phase_average(path+'\PhaseAveragePhi0_0002', 0, 2, 0, 900, 3)
data3 = phase_average(path+'\PhaseAveragePhi0_0003', 0, 3, 0, 900, 3)
data4 = phase_average(path+'\PhaseAveragePhi0_0004', 0, 4, 0, 900, 3)
#aveage data
data = np.zeros([NPZ*48, 9])
data = (data1 + data2 + data3 + data4)/4
os.chdir('d:\post')
outputfolder = 'post_result/'
#create output folder named 'post_result' 
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f = open( outputfolder + "phi0_newtest_fixed.plt",'w')
f.write("VARIABLES = X, Z, <U>/U*, <V>/U*, <W>/U*, <U'U'>/U*^2, <V'V'>/U*^2,\
         <W'W'>/U*^2, <-U'W'>/U*^2  \n")

np.savetxt(f, data)
f.close()

data5 = phase_average(path+'\PhaseAveragePhi1_0001', 1, 1, 0, 900, 1)
data6 = phase_average(path+'\PhaseAveragePhi1_0002', 1, 2, 0, 900, 1)
data7 = phase_average(path+'\PhaseAveragePhi1_0003', 1, 3, 0, 900, 1)
data8 = phase_average(path+'\PhaseAveragePhi1_0004', 1, 4, 0, 900, 1)
#aveage data
data0 = np.zeros([NPZ*48, 9])
data0 = (data5 + data6 + data7 + data8)/4

os.chdir('d:\post')
outputfolder = 'post_result/'
#create output folder named 'post_result' 
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f1 = open( outputfolder + "phi1_newtest_fixed.plt",'w')
f1.write("VARIABLES = X, Z, <U>/U*, <V>/U*, <W>/U*, <U'U'>/U*^2, <V'V'>/U*^2,\
         <W'W'>/U*^2, <-U'W'>/U*^2  \n")

np.savetxt(f1, data0)
f1.close()

data5 = phase_average(path+'\PhaseAveragePhi2_0001', 2, 1, 1, 31, 1)
data6 = phase_average(path+'\PhaseAveragePhi2_0002', 2, 2, 1, 31, 1)
data7 = phase_average(path+'\PhaseAveragePhi2_0003', 2, 3, 1, 31, 1)
data8 = phase_average(path+'\PhaseAveragePhi2_0004', 2, 4, 1, 31, 1)
#aveage data
data0 = np.zeros([NPZ*48, 9])
data0 = (data5 + data6 + data7 + data8)/4

os.chdir('d:\post')
outputfolder = 'post_result/'
#create output folder named 'post_result' 
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f1 = open( outputfolder + "phi2_fixed_new.plt",'w')
f1.write("VARIABLES = X, Z, <U>/U*, <V>/U*, <W>/U*, <U'U'>/U*^2, <V'V'>/U*^2,\
         <W'W'>/U*^2, <-U'W'>/U*^2  \n")

np.savetxt(f1, data0)
f1.close()

data9 = phase_average(path+'\PhaseAveragePhi3_0001', 3, 1, 1, 30, 1)
data10 = phase_average(path+'\PhaseAveragePhi3_0002', 3, 2, 1, 30, 1)
data11 = phase_average(path+'\PhaseAveragePhi3_0003', 3, 3, 1, 30, 1)
data12 = phase_average(path+'\PhaseAveragePhi3_0004', 3, 4, 1, 30, 1)
#aveage data
data0 = np.zeros([NPZ*48, 9])
data0 = (data9 + data10 + data11 + data12)/4

os.chdir('d:\post')
outputfolder = 'post_result/'
#create output folder named 'post_result' 
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f1 = open( outputfolder + "phi3_fixed_new.plt",'w')
f1.write("VARIABLES = X, Z, <U>/U*, <V>/U*, <W>/U*, <U'U'>/U*^2, <V'V'>/U*^2,\
         <W'W'>/U*^2, <-U'W'>/U*^2  \n")

np.savetxt(f1, data0)
f1.close()