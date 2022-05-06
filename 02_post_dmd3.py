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
modenum = 101 #how many modes will be used to reconstruct
iout1 = nx * (5 + 3) / 28.8 #turbine position : 5D, i output : 3D from turbine
jout1 = ny / 2
iout2 = nx * (5 + 5) / 28.8 #turbine position : 5D, i output : 3D from turbine
jout2 = ny / 2
iout1 = int(iout1)
jout1 = int(jout1)
iout2 = int(iout2)
jout2 = int(jout2)

casename = "114_2"
foldername = "./" + casename + "/POST_U_2D2_0004/"
#foldername = "./" + casename + "/POST_U_2D3_0001/"
outputfolder = "./" + casename + "_dmd/"

turbinefile = "./" + casename + "/" + "Turbine.inp"
turbinedata = np.genfromtxt(turbinefile, skip_header = 1)
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

f1 = open( outputfolder + "01.plt",'w')
f2 = open( outputfolder + "02_dmd_mode.plt",'w')
f3 = open( outputfolder + "03_reconstruction.plt",'w')
f4 = open( outputfolder + "04_comparison.plt",'w')

data1 = np.zeros((ibkb, ntm))
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

svdu0, svds, svdw_ct0 =np.linalg.svd(data1, full_matrices=False) #singular value decomposition # the column of svdu is pod modes.
svdu = np.zeros((ibkb,modenum))
svdw_ct = np.zeros((modenum,ntm))

svdu = svdu0[:,:modenum]
svdw_ct = svdw_ct0[:modenum,:]
svds_mat = np.diag(svds[:modenum])
svdu_c = np.transpose(svdu)  #conjugate transpose
svdu_ct = np.conjugate(svdu_c)
svdw_t = np.conjugate(svdw_ct) #conjugate transpose
svdw = np.transpose(svdw_t)
svds_mat_inv = np.linalg.inv(svds_mat)

svdst = np.dot(svdu_ct, np.dot(data2, np.dot(svdw,svds_mat_inv))) #reduced-order system matrix
evalue, evector = np.linalg.eig(svdst)

#phi = np.dot(data2, np.dot(svdw, np.dot(svds_mat_inv,evector))) #DMD modes
phi = np.dot(svdu, evector) #DMD modes #Scmid 2010 JFM
logeval = np.log(evalue)
revalue = np.real(evalue)
ievalue = np.imag(evalue)
omega = np.imag(logeval)/dt
sigma = np.real(logeval)/dt

scaling_coeff_dmd = np.zeros((modenum), dtype = complex)
scaled_phi = np.zeros((ibkb,modenum), dtype = complex)
vandermonde_mat = np.zeros((modenum, ntm), dtype = complex)
recon_dmd = np.zeros((ibkb,ntm), dtype = complex) #reconstruction flow field

scaling_coeff_dmd = np.dot(np.linalg.pinv(phi[:,:modenum]),data1[:,0])
scaling_coeff_mat = np.diag(scaling_coeff_dmd[:modenum])
for it in range(ntm):
	vandermonde_mat[:,it] = np.power(evalue[:],it, dtype = complex)
recon_dmd = np.dot(phi, np.dot(scaling_coeff_mat, vandermonde_mat))
for it in range(modenum):
	scaled_phi[:,it] = scaling_coeff_dmd[it]*phi[:,it]
data3 = np.zeros((modenum,7))
for i in range(modenum):
	data3[i,0] = i
	data3[i,1] = svds[i]
	data3[i,2] = np.absolute(scaling_coeff_dmd[i])
	data3[i,3] = revalue[i]
	data3[i,4] = ievalue[i]
	data3[i,5] = sigma[i]
	data3[i,6] = omega[i]

f1.write("VARIABLES = IT, singular_value,abs(scaling_coef), ev_r, ev_i, ln(ev_r)(sigma), ln(ev_i)(omega) \n")
np.savetxt(f1,data3)
f2.write("VARIABLES = X, Y, DMD_Mode, Scaled_DMD_Mode \n")
f3.write("VARIABLES = X, Y, Recon_dmd, Original \n")
f4.write("VARIABLES = time, ReconDmd_pos1, Original_pos1, ReconDmd_pos2, Original_pos2 \n")
data4 = np.zeros((ibkb,4))
data5 = np.zeros((ibkb,4))
data6 = np.zeros((ntm,5))

for it in range(modenum):
	f2.write("ZONE T = '" + str(it) + "' I = " + str(ib2-ib1+1) + " J = " + str(kb2-kb1+1) + "\n" )
	d = 0
	for k in range(kb1,kb2+1):
		for i in range(ib1,ib2+1):
			j = (k-kb1)*(ib2-ib1+1)+(i-ib1)
			data4[j,0] = x1[i]
			data4[j,1] = z1[k]
			data4[j,2] = phi[d,it]
			data4[j,3] = scaled_phi[d,it]
			d = d + 1
	np.savetxt(f2,data4)

for it in range(ntm): #originally range(ntm)
	if it < 50:
		f3.write("ZONE T = '" + str(it) + "' I = " + str(ib2-ib1+1) + " J = " + str(kb2-kb1+1) + "\n" )
	d = 0
	for k in range(kb1,kb2+1):
		for i in range(ib1,ib2+1):
			j = (k-kb1)*(ib2-ib1+1)+(i-ib1)
			data5[j,0] = x1[i]
			data5[j,1] = z1[k]
			data5[j,2] = recon_dmd[d,it]
			data5[j,3] = data1[d,it]
			if i == iout1 and k == jout1:
				data6[it,0] = it * dt
				data6[it,1] = recon_dmd[d,it]
				data6[it,2] = data1[d,it]
			if i == iout2 and k == jout2:
				data6[it,3] = recon_dmd[d,it]
				data6[it,4] = data1[d,it]
			d = d + 1

	if it < 50:
		np.savetxt(f3,data5)
np.savetxt(f4,data6)

f1.close()
f2.close()
f3.close()
f4.close()

#	UM = np.zeros((nx,ny))
#	for iit in range(it - ctsi * noutc , it + ctsi * noutc + 1, noutc):
#		filename = foldername + "POST_U_2D2_" + "{:010d}".format(iit) + "_0004.DAT"
#		data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
#		data0 = data0.reshape([ny,nx,nvar])
#		U = data0[:,:,3].reshape([ny,nx]).transpose()
#		UM = UM + U / (ctsi * 2 + 1)
	
#	if it < tis + 300:
#		nline3 = nx*ny
#		data3 = np.zeros((nline3,3))
#		for j in range(ny):
#			for i in range(nx):
#				k = j*nx + i
#				data3[k,0] = x1[i]
#				data3[k,1] = y1[j]
#				data3[k,2] = UM[i,j]

#		f2 = open( outputfolder + "POST_UM_2D2_" + "{:010d}".format(it) + ".DAT" , 'w')
#		f2.write("VARIABLES = X, Y, UM \n")
#		f2.write("Zone T =" + "'" + str(it) + "'" + " I = " + str(nx) + " J = " + str(ny) + "\n")
#		np.savetxt(f2,data3)
