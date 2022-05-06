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

def identify_wake_core(_ixran, _iyran, _h, _x1, _y1, width=20, _method=2, _hw=6):
	_ix = []
	_iy = []
	_x = []
	_y = []
	#print("ixran="+str(_ixran)+", iyran="+str(_iyran)+", iyran[0]="+str(_iyran[0]))
	
	_tempwidth = (_iyran[1] - _iyran[0])/3
	_tempstart = _iyran[0] + (_iyran[1] - _iyran[0])/2 - _tempwidth/2	
	_tempend = _tempstart + _tempwidth
	
	_utmp = np.zeros(_h.shape[1])
	
	for _i in range(_ixran[0], _ixran[1]+1):
		#_ihmin = np.argmin(_h[_i, _iyran[0]:(_iyran[1]+1)])
		#print("ihmin="+str(_ihmin))
		#_itemp = _iyran[0] + _ihmin
	
		if _method==1:
			# method 1: simply find the global minimum
			_ihmin = np.argmin(_h[_i, _tempstart:_tempend]) # ... .. ... index . .....
		elif _method==2:
			# method 2: firstly get a band average value, then find the location of the minimum
			for _j in range(_iyran[0]+_hw, _iyran[1]-_hw):
				_utmp[_j] = np.average(_h[_i, (_j-_hw):(_j+_hw)]) # y.... smoothing . . ... minimum index. ....
			_ihmin = np.argmin(_utmp[_tempstart:_tempend])
			
		## I'd like to set the hub as starting point of wake centerline
				
		
		## following code is to improve the continuity of search (only search for the next point in a range of last point.)
		## but it turns out to be easily limited to bottom or top edge of the search frame		
		_itemp = _tempstart + _ihmin
		_tempstart = _itemp - width/2
		
		_tempwidth = width
		
		_tempend = _itemp + width/2
		if _tempstart < (_iyran[0]+_hw):
			_tempstart = _iyran[0] + _hw	
		if _tempend > (_iyran[1]-_hw):
			_tempend = _iyran[1] - _hw
		
		_ix.append(_i)
		_iy.append(_itemp)
		_x.append(_x1[_i])
		_y.append(_y1[_itemp])
		#print(str(_i)+" "+str(_itemp)+" "+str(_x2[_i])+" "+str(_y2[_itemp]))
	return _ix, _iy, _x, _y

def smooth(x,window_len=7,window='flat'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]  #np.r_ : .. ... (row ...)
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')    #np.hanning(11) ... .... .. . .... w. ..... --- hanning window : pyhon .. window.

    y=np.convolve(w/w.sum(),s,mode='valid')
    #print(len(y))
    #print(len(w))
    #return y
    
    _start = len(w)/2
    return y[(_start):(_start+len(x))]

r = 0.075
D = 2.0*r
nx = 384
ny = 96
nz = 1
nvar = 6
nturbine = 1
tiss = 10000
tiee = 30000
nwriteperperiod = 100
coherencetimescale = 0.6 #0.6 * T
noutc = 2
#coherence time scale index --- number of averaging in one side
ctsi = nwriteperperiod * coherencetimescale / noutc / 2
ctsi = int(ctsi)

casename = "114_2"

foldername = "./" + casename + "/POST_U_2D2_0004/"
outputfolder = "./" + casename + "_wake/"

tis = tiss + ctsi * noutc
tie = tiee - ctsi * noutc + 1
tii = noutc

y_ts = -2.0/3.0*r

yfftm2 = []
u_ts2 = []

turbinefile = "./" + casename + "/" + "Turbine.inp"
turbinedata = np.genfromtxt(turbinefile, skip_header = 1)
turbineloc = turbinedata[3:5] #single turbine
#nturbine = len(turbineloc[:,0]) #multiple turbines

for iturb in range(nturbine):
	xtmp, ytmp = get_coord_range(turbineloc[:],D)

allmins = []
for iturb in range(nturbine):
	allmins.append(list())

f = open( outputfolder + "wake_spectrum.dat",'w')
f1 = open( outputfolder + "wake_line.dat",'w')

for it in range(tis, tie, tii):
	filename = foldername + "POST_U_2D2_" + "{:010d}".format(it) + "_0004.DAT"
	print(filename)
	data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
	data0 = data0.reshape([ny,nx,nvar])
	X = data0[:,:,0].reshape([ny,nx]).transpose()
	Y = data0[:,:,1].reshape([ny,nx]).transpose()
#	U = data0[:,:,3].reshape([ny,nx]).transpose()
	x1 = X[:,0]
	y1 = Y[0,:]
	
	UM = np.zeros((nx,ny))
	for iit in range(it - ctsi * noutc , it + ctsi * noutc + 1, noutc):
		filename = foldername + "POST_U_2D2_" + "{:010d}".format(iit) + "_0004.DAT"
		data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
		data0 = data0.reshape([ny,nx,nvar])
		U = data0[:,:,3].reshape([ny,nx]).transpose()
		UM = UM + U / (ctsi * 2 + 1)
	
	if it < tis + 100:
		nline3 = nx*ny
		data3 = np.zeros((nline3,3))
		for j in range(ny):
			for i in range(nx):
				k = j*nx + i
				data3[k,0] = x1[i]
				data3[k,1] = y1[j]
				data3[k,2] = UM[i,j]

		f2 = open( outputfolder + "POST_UM_2D2_" + "{:010d}".format(it) + ".DAT" , 'w')
		f2.write("VARIABLES = X, Y, UM \n")
		f2.write("Zone T =" + "'" + str(it) + "'" + " I = " + str(nx) + " J = " + str(ny) + "\n")
		np.savetxt(f2,data3)

	xran = np.zeros((nturbine,2))
	yran = np.zeros((nturbine,2))
	ixran = np.zeros((nturbine,2), dtype = 'int')
	iyran = np.zeros((nturbine,2), dtype = 'int')
	xwake = []
	ywake = []
	ywpo = []
	ywp = []
	yfft = []
	yfftm = []
	u_ts = []

	for iturb in range(nturbine):
#		xran[iturb,:], yran[iturb,:] = get_coord_range(turbineloc[iturb,:], D)
#		ixran[iturb,:], iyran[iturb,:] = get_index_range(x2, y2, xran[iturb,:], yran[iturb,:])
		xran[iturb,:], yran[iturb,:] = get_coord_range(turbineloc[:], D)
		ixran[iturb,:], iyran[iturb,:] = get_index_range(x1, y1, xran[iturb,:], yran[iturb,:])
		y2_ts = (yran[iturb,0] + yran[iturb,1])/2.0 + y_ts

		ixtemp,iytemp,xtemp,ytemp = identify_wake_core(ixran[iturb,:],iyran[iturb,:],UM,x1,y1)
		xtemp = np.array(xtemp)
		ytemp = np.array(ytemp)
		xwake.append(xtemp)
		ywake.append(ytemp)
		ywpo.append(ytemp)

		np_tmp = np.empty((len(ytemp), 2))
		np_tmp[:,0] = xtemp
		np_tmp[:,1] = ytemp
		allmins[iturb].append(np_tmp)
			
		ywp = smooth(ytemp)
		#print("after smooth: len=" + str(len(ytemp)))
		#ywp.append(ytemp)
		#print("iturb="+str(iturb))
		#print(np.array(ixtemp))
		#print(np.array(iytemp))
		#print(np.vstack((np.array(ixtemp), np.array(iytemp))))

		p_temp = np.abs(np.fft.fft(ytemp - np.average(ytemp)))
		yfft.append(p_temp)
		if iturb==0:
			yfftm = p_temp
		else:
			yfftm = yfftm + p_temp / nturbine

		xfft = 2.0*np.pi/(xran[iturb,1]-xran[iturb,0]) * np.arange(len(yfftm))

		nline2 = len(xtemp)
		data2 = np.zeros((nline2,3))
		for i in range(nline2):
			data2[i,0] = xtemp[i]
			data2[i,1] = ytemp[i]
			data2[i,2] = ywp[i]
		
	if it==tis:
		yfftm2 = yfftm
	else:
		yfftm2 = yfftm2 + yfftm / len(range(tis, tie, tii))
	
	nline = len(xfft)
	data1 = np.zeros((nline,3))
	for i in range(nline):
		data1[i,0] = xfft[i]
		data1[i,1] = yfftm[i]
		data1[i,2] = yfftm2[i]

	f.write("Zone T =" + "'" + str(it) + "'\n")
	np.savetxt(f,data1)
	f1.write("Zone T =" + "'" + str(it) + "'\n")
	np.savetxt(f1,data2)
	#np.savetxt( casename + "_wake_spectrum.dat", data1)

f.close()
f1.close()
