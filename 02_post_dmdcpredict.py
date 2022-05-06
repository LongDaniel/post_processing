import sys
import numpy as np
import tecplot_io as tec

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

def get_coord_range (_loc, _D):
	_xc = _loc[0]
	_yc = _loc[1]
	_xran = [_xc-3.0*_D, _xc+9.0*_D]
	_yran = [_yc-2.0*_D, _yc+2.0*_D]
	return _xran, _yran

def get_index_range (_x, _y, _xran, _yran):
#	print "range for x="+str(_x)+" is "+str(_xran)
#	print "range for y="+str(_y)+" is "+str(_yran)
	_ixmin, _ixmax = find_index(_x, _xran)
	_iymin, _iymax = find_index(_y, _yran)
	_ixran = [_ixmin, _ixmax]
	_iyran = [_iymin, _iymax]
	return _ixran, _iyran

def get_coord_range_x (_loc, _D):
	_xc = _loc[0]
	_xran = [_xc-3.0*_D, _xc+9.0*_D]
	return _xran

def get_index_range_x (_x, _xran):
#	print "range for x="+str(_x)+" is "+str(_xran)
#	print "range for y="+str(_y)+" is "+str(_yran)
	_ixmin, _ixmax = find_index(_x, _xran)
	_ixran = [_ixmin, _ixmax]
	return _ixran

r = 0.075
D = 2.0*r
nx = 384
ny = 96
nz = 1
#nz = 65
nvar = 6
nturbine = 1
tis = 10020
tie = 20000 + 1
tii = 20
dt = 0.0013662 * tii
nt = (tie -1 - tis) / tii + 1
ntm = nt - 1
modenum = 150  #truncation value of data 1 SVD
modenum2 = 149 #truncation value of data 2 SVD
ipredict = 1
iout1 = nx * (5 + 3) / 28.8 #turbine position : 5D, i output : 3D from turbine
jout1 = ny / 2
iout2 = nx * (5 + 5) / 28.8 #turbine position : 5D, i output : 3D from turbine
jout2 = ny / 2
iout1 = int(iout1)
jout1 = int(jout1)
iout2 = int(iout2)
jout2 = int(jout2)

casename_system = "114_4"
casename = "114_4"
foldername = "./" + casename + "/POST_U_2D2_0004/"
#foldername = "./" + casename + "/POST_U_2D3_0001/"
outputfolder = "./" + casename + "_dmd_predict/"
systemfolder = "./" + casename_system + "_dmd/"


turbinefile = "./" + casename + "/" + "Turbine.inp"
turbinedata = np.genfromtxt(turbinefile, skip_header = 1)
bladepitchcontrol = "./" + casename + "/" +"Turbine_AL_001.dat"
bladepitchdata = np.genfromtxt(bladepitchcontrol, skip_header = 1)
nbladepitchdata = len(bladepitchdata[:,0])

bpa = np.zeros(ntm)
d = 0
for i in range(tis, tie-1, tii):
	bpa[d] = bladepitchdata[i,8]
	d = d + 1

turbineloc = turbinedata[3:5] #single turbine
#nturbine = len(turbineloc[:,0]) #multiple turbines

for iturb in range(nturbine):
	xtmp, ytmp = get_coord_range(turbineloc[:],D)

filename = foldername + "POST_U_2D2_" + "{:010d}".format(tis) + "_0004.DAT"
#filename = foldername + "POST_U_2D3_" + "{:010d}".format(tis) + "_0001.DAT"
data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
#x-y direc
data0 = data0.reshape([ny,nx,nvar])
X = data0[:,:,0].reshape([ny,nx]).transpose()
Z = data0[:,:,1].reshape([ny,nx]).transpose() #be careful the index of data0 [:,:,1 - Y value, :,:,2 - Z value]
#x-z direc
#data0 = data0.reshape([nz,nx,nvar])
#X = data0[:,:,0].reshape([nz,nx]).transpose()
#Z = data0[:,:,2].reshape([nz,nx]).transpose()
x1 = X[:,0]
z1 = Z[0,:]

xran = np.zeros((nturbine,2))
yran = np.zeros((nturbine,2))
ixran = np.zeros((nturbine,2), dtype = 'int')
iyran = np.zeros((nturbine,2), dtype = 'int')
for iturb in range(nturbine):
#x-y direc
	xran[iturb,:], yran[iturb,:] = get_coord_range(turbineloc[:], D)
	ixran[iturb,:], iyran[iturb,:] = get_index_range(x1, z1, xran[iturb,:], yran[iturb,:])
#x-z dir
#	xran[iturb,:] = get_coord_range_x(turbineloc[:], D)
#	ixran[iturb,:] = get_index_range_x(x1, xran[iturb,:])

ib1 = ixran[0,0]
ib2 = ixran[0,1]
#x-y direc
kb1 = iyran[0,0]
kb2 = iyran[0,1]
#x-z direc
#kb1 = 0
#kb2 = nz - 3
ibkb = (ib2 - ib1 + 1) * (kb2 - kb1 + 1)

data1 = np.zeros((ibkb + 1, ntm))
data2 = np.zeros((ibkb, ntm))
for it in range(tis, tie, tii):
	it1 = (it - tis) / tii 
	filename = foldername + "POST_U_2D2_" + "{:010d}".format(it) + "_0004.DAT"
#	filename = foldername + "POST_U_2D3_" + "{:010d}".format(it) + "_0001.DAT"
	print(filename)
	data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
#x-y direc
	data0 = data0.reshape([ny,nx,nvar])
	U = data0[:,:,3].reshape([ny,nx]).transpose()
#x-z direc
#	data0 = data0.reshape([nz,nx,nvar])
#	U = data0[:,:,3].reshape([nz,nx]).transpose() #x-z direc
	dd = 0
	for k in range (kb1, kb2 + 1):
		for i in range(ib1, ib2 + 1):
			if it1 <= ntm - 1:
				data1[dd, it1] = U[i,k] #data set from 1 to n-1
			if it1 >= 1:
				data2[dd, it1-1] = U[i,k] #data set at n
			dd = dd + 1
	if it1 == ntm:
		data1[dd , :] = bpa[:]


#reduced-order system matrix read & output########################################
f3 = open( outputfolder + "03_reconstruction.plt",'w')
f4 = open( outputfolder + "04_comparison.plt",'w')
f5 = open( systemfolder + "00_svdu_ct_.dat",'r')
f6 = open( systemfolder + "00_Atil.dat",'r')
f7 = open( systemfolder + "00_Btil.dat",'r')
svdu_ct_ = np.zeros((modenum2, ibkb))
Atil = np.zeros((modenum2, modenum2))
Btil = np.zeros(modenum2)
svdu_ct_read = np.loadtxt(f5)
Atil = np.loadtxt(f6)
Btil = np.loadtxt(f7)
f5.close()
f6.close()
f7.close()

recon_dmdc = np.zeros(ntm, dtype = complex)
recon_dmdc_wocon = np.zeros(ntm, dtype = complex)
recon_dmdc_full = np.zeros((ibkb,ntm), dtype = complex)
recon_dmdc_full_wocon = np.zeros((ibkb,ntm), dtype = complex)

for it in range(ntm):
		if it == 0:
			recon_dmdc_full[:,0] = data2[:,0]
			recon_dmdc_full_wocon[:,0] = data2[:,0]
			recon_dmdc = np.dot(svdu_ct_, recon_dmdc_full[:,0])
			recon_dmdc_wocon = np.dot(svdu_ct_, recon_dmdc_full_wocon[:,0]) #x(k+1)til = svdu_ct_ * x(k+1) : coordinate transformation based on POD modes
		else:
			recon_dmdc = np.dot(Atil, recon_dmdc) + Btil[:]*bpa[it]  #x(k+1)til = Atil * x(k)til + Btil * u(k)
			recon_dmdc_wocon = np.dot(Atil, recon_dmdc_wocon)
			recon_dmdc_full[:,it] = np.dot(svdu_, recon_dmdc)
			recon_dmdc_full_wocon[:,it] = np.dot(svdu_, recon_dmdc_wocon)


f3.write("VARIABLES = X, Y, Recon_dmd, Recon_dmdc, Original \n")
f4.write("VARIABLES = time, ReconDmd1,ReconDmdC1_cx,ReconDmdC1_c,Original1,ReconDmd2,ReconDmdC2_cx,ReconDmdC2_c,Original2,input(bpa) \n")

data5 = np.zeros((ibkb,5))     #f3
data6 = np.zeros((ntm,10))     #f4

#f3&f4##########################################
for it in range(ntm): #originally range(ntm)
	if it < 50:
		f3.write("ZONE T = '" + str(it) + "' I = " + str(ib2-ib1+1) + " J = " + str(kb2-kb1+1) + "\n" )
	d = 0
	for k in range(kb1,kb2+1):
		for i in range(ib1,ib2+1):
			j = (k-kb1)*(ib2-ib1+1)+(i-ib1)
			data5[j,0] = x1[i]
			data5[j,1] = z1[k]
			data5[j,2] = recon_dmd[d,it]          #data reconstruction only via a reduced-order approximation of the system matrix A
			data5[j,3] = recon_dmdc_full[d,it]
			data5[j,4] = data2[d,it]
			if i == iout1 and k == jout1:
				data6[it,0] = it * dt
				data6[it,1] = recon_dmd[d,it] #data reconstruction only via a reduced-order approximation of the system matrix A
				data6[it,2] = recon_dmdc_full_wocon[d,it]
				data6[it,3] = recon_dmdc_full[d,it]
				data6[it,4] = data2[d,it]
			if i == iout2 and k == jout2:
				data6[it,5] = recon_dmd[d,it] #data reconstruction only via a reduced-order approximation of the system matrix A
				data6[it,6] = recon_dmdc_full_wocon[d,it]
				data6[it,7] = recon_dmdc_full[d,it]
				data6[it,8] = data2[d,it]
				data6[it,9] = bpa[it]
			d = d + 1

	if it < 50:
		np.savetxt(f3,data5)
np.savetxt(f4,data6)

f3.close()
f4.close()

