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
	_xran = [_xc+0.1*_D, _xc+8.0*_D]
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
ny = 1
nz = 65
nvar = 6
nturbine = 1
tis = 10000
tie = 10290 + 1
tii = 10
dt = 0.0013662*tii
nt = (tie -1 - tis) / tii + 1
ntm = nt - 1
modenum = 11 #how many modes will be used to reconstruc

casename = "114_2"
foldername = "./" + casename + "/POST_U_2D3_0001/"
outputfolder = "./" + casename + "_dmd/"

turbinefile = "./" + casename + "/" + "Turbine.inp"
turbinedata = np.genfromtxt(turbinefile, skip_header = 1)
turbineloc = turbinedata[3:5] #single turbine
#nturbine = len(turbineloc[:,0]) #multiple turbines

for iturb in range(nturbine):
	xtmp, ytmp = get_coord_range(turbineloc[:],D)

filename = foldername + "POST_U_2D3_" + "{:010d}".format(tis) + "_0001.DAT"
data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
data0 = data0.reshape([nz,nx,nvar])
X = data0[:,:,0].reshape([nz,nx]).transpose()
#Y = data0[:,:,1].reshape([ny,nx]).transpose()
Z = data0[:,:,2].reshape([nz,nx]).transpose()
x1 = X[:,0]
z1 = Z[0,:]

xran = np.zeros((nturbine,2))
yran = np.zeros((nturbine,2))
ixran = np.zeros((nturbine,2), dtype = 'int')
iyran = np.zeros((nturbine,2), dtype = 'int')
for iturb in range(nturbine):
#	xran[iturb,:], yran[iturb,:] = get_coord_range(turbineloc[:], D)
#	ixran[iturb,:], iyran[iturb,:] = get_index_range(x1, y1, xran[iturb,:], yran[iturb,:])
	xran[iturb,:] = get_coord_range_x(turbineloc[:], D)
	ixran[iturb,:] = get_index_range_x(x1, xran[iturb,:])

ib1 = ixran[0,0]
ib2 = ixran[0,1]
kb1 = 0
kb2 = nz - 3
ibkb = (ib2 - ib1 + 1) * (kb2 - kb1 + 1)

f1 = open( outputfolder + "01.dat",'w')
f2 = open( outputfolder + "02.dat",'w')

data1 = np.zeros((ibkb, ntm))
data2 = np.zeros((ibkb, ntm))
for it in range(tis, tie, tii):
	it1 = (it - tis) / tii 
	filename = foldername + "POST_U_2D3_" + "{:010d}".format(it) + "_0001.DAT"
	print(filename)
	data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
	data0 = data0.reshape([nz,nx,nvar])
#	X = data0[:,:,0].reshape([nz,nx]).transpose()
#	Y = data0[:,:,1].reshape([ny,nx]).transpose()
#	Z = data0[:,:,2].reshape([nz,nx]).transpose()
	U = data0[:,:,3].reshape([nz,nx]).transpose()
#	x1 = X[:,0]
#	z1 = Z[0,:]
	dd = 0
	for k in range (kb1, kb2 + 1):
		for i in range(ib1, ib2 + 1):
			if it1 <= ntm - 1:
				data1[dd, it1] = U[i,k] #data set from 1 to n-1
			if it1 >= 1:
				data2[dd, it1-1]      = U[i,k] #data set at n
			dd = dd + 1

svdu, svds, svdw_ct =np.linalg.svd(data1, full_matrices=False) #singular value decomposition # the column of svdu is pod modes.
svds_mat = np.diag(svds)
svdu_c = np.transpose(svdu)  #conjugate transpose
svdu_ct = np.conjugate(svdu_c)
svdw_t = np.conjugate(svdw_ct) #conjugate transpose
svdw = np.transpose(svdw_t)
svds_mat_inv = np.linalg.inv(svds_mat)
svdst = np.dot(svdu_ct, np.dot(data2, np.dot(svdw,svds_mat_inv)))
evalue, evector = np.linalg.eig(svdst)
logeval = np.log(evalue)
revalue = np.real(evalue)
ievalue = np.imag(evalue)
omega = np.imag(logeval)/dt
sigma = np.real(logeval)/dt
#omega = np.angle(eigen value)/dt
#sigma = np.log(np.absolute(eigen value))/dt
phi = np.dot(data2, np.dot(svdw, np.dot(svds_mat_inv,evector)))
normphi = np.zeros(ntm)
normsvdu = np.zeros(ntm)
scaling_coeff_dmd = np.zeros((modenum), dtype = complex)
scaling_coeff_dmd = np.dot(np.linalg.pinv(phi[:,:modenum]),data1[:,0])
proj_coeff_pod = np.zeros((ntm,ntm))
recon_dmd = np.zeros((ibkb,ntm), dtype = complex)
recon_pod = np.zeros((ibkb,ntm))
for it in range(ntm):
	normphi[it] = np.linalg.norm(phi[:,it])
	normsvdu[it] = np.linalg.norm(svdu[:,it])

for it in range(ntm):
	for it1 in range(modenum):   #How many modes you want to use to reconstruct the original flow field?
		proj_coeff_pod[it, it1] = np.inner(svdu[:,it1],data1[:,it])/(normsvdu[it1]**2.)
		recon_dmd[:,it] = recon_dmd[:,it] + scaling_coeff_dmd[it1]*np.power(evalue[it1],it)*phi[:,it1]
		recon_pod[:,it] = recon_pod[:,it] + proj_coeff_pod[it, it1]*svdu[:,it1]

data3 = np.zeros((ntm,7))
for i in range(ntm):
	data3[i,0] = i
	data3[i,1] = revalue[i]
	data3[i,2] = ievalue[i]
	data3[i,3] = sigma[i]
	data3[i,4] = omega[i]
	data3[i,5] = normphi[i]
	data3[i,6] = svds[i]

f1.write("VARIABLES = IT, Real_eigen, Imag_eigen, LogReal_eigen(sigma), LogImag_eigen(omega), norm, sigular_value \n")
np.savetxt(f1,data3)
f2.write("VARIABLES = X, Y, DMD_Mode, POD_Mode, Recon_dmd, Recon_pod, Original \n")
data4 = np.zeros((ibkb,7))
for it in range(ntm):
	f2.write("ZONE T = '" + str(it) + "' I = " + str(ib2-ib1+1) + " J = " + str(kb2-kb1+1) + "\n" )
	d = 0
	for k in range(kb1,kb2+1):
		for i in range(ib1,ib2+1):
			j = (k-kb1)*(ib2-ib1+1)+(i-ib1)
			data4[j,0] = x1[i]
			data4[j,1] = z1[k]
			data4[j,2] = phi[d,it]
			data4[j,3] = svdu[d,it]
			data4[j,4] = recon_dmd[d,it]
			data4[j,5] = recon_pod[d,it]
			data4[j,6] = data1[d,it]
			d = d + 1
	np.savetxt(f2,data4)

f1.close()
f2.close()



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
