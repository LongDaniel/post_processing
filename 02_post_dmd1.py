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
	_xran = [_xc+0.2*_D, _xc+23.0*_D]
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
	_xran = [_xc+0.05*_D, _xc+23.0*_D]
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
tie = 11500 + 1
tii = 30
dt = 0.0013662
nt = (tie -1 - tis) / tii + 1
ntm = nt - 1

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
f3 = open( outputfolder + "03.dat",'w')

data1 = np.zeros((ibkb, ntm))
data2 = np.zeros(ibkb)
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
			else:
				data2[dd]      = U[i,k] #data set at n
			dd = dd + 1

data1T = np.transpose(data1)
aa = np.dot(data1T, data1) #ntm*ntm /left inverse of non-square matrix A is (A^T A)^-1 A^T
aai = np.linalg.inv(aa)    #ntm*ntm
ali = np.dot(aai,data1T)   #ntm*ibkb /left inverse of a (data1)
c = np.dot(ali,data2)      #ntm*1    /v^n = v^(1 to n-1) c
com = np.zeros((ntm,ntm))
for ii in range(ntm-1):
	for jj in range(ntm-1):
		if ii == jj:
			com[ii+1,jj] = 1 #ntm * ibkb
com[:,ntm-1] = c[:]
evalue, evector = np.linalg.eig(com)
logeval = np.log(evalue)
revalue = np.real(evalue)
ievalue = np.imag(evalue)
omega = np.imag(logeval)/dt
sigma = np.real(logeval)/dt
#omega = np.angle(evalue)/dt
#sigma = np.log(np.absolute(evalue))/dt
phi = np.dot(data1,evector) #ibkb*ntm
normphi = np.zeros((ntm))
for it in range(ntm):
	normphi[it] = np.linalg.norm(phi[:,it])

data3 = np.zeros((ntm,6))
for i in range(ntm):
	data3[i,0] = i
	data3[i,1] = revalue[i]
	data3[i,2] = ievalue[i]
	data3[i,3] = sigma[i]
	data3[i,4] = omega[i]
	data3[i,5] = normphi[i]

f1.write("VARIABLES = IT, Real_eigen, Imag_eigen, LogReal_eigen(sigma), LogImag_eigen(omega), norm \n")
np.savetxt(f1,data3)
f2.write("VARIABLES = X, Y, Mode, Original \n")
np.savetxt(f3,evalue)
np.savetxt(f3,evalue.shape)
np.savetxt(f3,evector)
np.savetxt(f3,evector.shape)
np.savetxt(f3,phi)
np.savetxt(f3,phi.shape)
data4 = np.zeros((ibkb,4))
for it in range(ntm):
	f2.write("ZONE T = '" + str(it) + "' I = " + str(ib2-ib1+1) + " J = " + str(kb2-kb1+1) + "\n" )
	d = 0
	for k in range(kb1,kb2+1):
		for i in range(ib1,ib2+1):
			j = (k-kb1)*(ib2-ib1+1)+(i-ib1)
			data4[j,0] = x1[i]
			data4[j,1] = z1[k]
			data4[j,2] = phi[d,it]
			data4[j,3] = data1[d,it]
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

