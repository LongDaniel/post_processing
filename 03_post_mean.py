import sys
import h5py
import numpy as np
import tecplot_io as tec
import os

def func2(x, a, b):
  return a*x+b

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

#kappa = 0.4
#nu = 1.511e-5
#PEX = 1.45444104333
#PEY = 8.72664625997
#hbar = 0.46
#uinfty = 2.54390548295
dt = 0.68543297937
#Rotational angular period
T_turb = 42.84
U_star = 0.356
H_hub = 70
#Mean finite velocity
U = 11.5258407161

print("Current working directory: {0}".format(os.getcwd()) )

#working directory path
path = 'd:\post\Fixed'
os.chdir(path)

outputfolder = 'post_result/'
#create output folder named 'post_result' 
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
print("Current working directory: {0}".format(os.getcwd()))

tis = 150
tie = 45000
tii = 150
nti = int((tie - tis) / tii + 1)

f3 = open( outputfolder + "01_um,vm,wm,uum,vvm,wwm.plt",'w')
f4 = open( outputfolder + "01_utm,vtm,wtm,uutm,vvtm,wwtm.plt",'w')
f5 = open( outputfolder + "03_mean_flux.plt",'w')
f3.write("VARIABLES = z/H_hub, um/u*, vm, wm, uum, vvm, wwm, uvm, uwm, vwm  \n")
f4.write("VARIABLES = z/H_hub, utm/u*, vtm, wtm, uutm, vvtm, wwtm, uvtm, uwtm, vwtm  \n")
f5.write("VARIABLES = t/T, mean_flux \n")

for it in range(nti):
    ti = tis + tii * it
    time = ti * dt
    f3.write("ZONE T = '" + str(time) + "'" + "\n" )
    fname = 'DAT_{:010d}.h5'.format(ti)
#    print("Reading file "+ fname)
    f = h5py.File(fname, "r")

    #print("Keys: %s" % f.keys())
    ## Old version Keys: [u'dz', u'dzw', u'eta', u'eta0', u'hh', u'pp', u'u', u'v', u'w', u'z', u'zw', u'zz']
    ## New version Keys: [u'eta', u'hh', u'pp', u'u', u'v', u'w', u'z']
    
    zz = np.array(f["z"][:,0,0]).copy()
    u = f["u"]
    v = f["v"]
    w2 = f["w"]
    w = np.array(w2).copy()

#   print(u.shape)

    NPX = u.shape[2]
    NPY = u.shape[1]
    NPZ = u.shape[0]
    
    if it==0:
      data3 = np.zeros((NPZ, 10))
      data4 = np.zeros((NPZ, 10))
      data5 = np.zeros((nti, 2))
      u_m_all = np.zeros((NPZ,3))
      uu_m_all = np.zeros((NPZ,6))

    w[0, :, :] = w2[0, :, :]
    for k in range(1,NPZ):
      w[k, :, :] = 0.5*(w2[k-1, :, :] + w2[k, :, :])

    u_m = np.zeros((NPZ,3))
    uu_m1 = np.zeros((NPZ,6))
    uu_m2 = np.zeros((NPZ,6))

    for k in range(NPZ):
      u_m[k, 0] = np.average(u[k,:,:])
      u_m[k, 1] = np.average(v[k,:,:])
      u_m[k, 2] = np.average(w[k,:,:])
      uu_m1[k, 0] = np.average(u[k,:,:]**2)
      uu_m1[k, 1] = np.average(v[k,:,:]**2)
      uu_m1[k, 2] = np.average(w[k,:,:]**2)
      uu_m1[k, 3] = np.average(u[k,:,:]*v[k,:,:])
      uu_m1[k, 4] = np.average(u[k,:,:]*w[k,:,:])
      uu_m1[k, 5] = np.average(v[k,:,:]*w[k,:,:])
      uu_m2[k, 0:3] = uu_m1[k, 0:3] - u_m[k,0:3]**2
      uu_m2[k, 3] = uu_m1[k, 3] - u_m[k, 0] * u_m[k, 1]
      uu_m2[k, 4] = uu_m1[k, 4] - u_m[k, 0] * u_m[k, 2]
      uu_m2[k, 5] = uu_m1[k, 5] - u_m[k, 1] * u_m[k, 2]
      data3[k, 0] = zz[k] / H_hub
      data3[k, 1:4] = (u_m[k,0:3] * U) / U_star
      data3[k, 4:10] = uu_m2[k,0:6]
    np.savetxt(f3,data3)

    u_m_all = u_m_all + u_m
    uu_m_all = uu_m_all + uu_m2
    f.close()

    mean_flux = np.trapz(u_m[:,0], zz[:])
    data5[it, 0] = time / T_turb
    data5[it, 1] = mean_flux
    print(str(ti) +"\t"+ str(mean_flux))

u_m_all = u_m_all / nti
uu_m_all = uu_m_all / nti

for k in range(NPZ):
    data4[k,0] = zz[k]
    data4[k,1:4] = u_m_all[k,0:3]
    data4[k,4:10] = uu_m_all[k,0:6]
np.savetxt(f4,data4)
np.savetxt(f5,data5)
f3.close()
f4.close()
f5.close()
