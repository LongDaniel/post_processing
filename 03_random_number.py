import sys
import numpy as np
import random
import tecplot_io as tec

dur = 100
tot = 30000
ran = tot/dur

casename = "114_2"
outputfolder = "./" + casename

#f1 = open( outputfolder + "bladepitching.inp",'w')
f1 = open("./" + "bladepitching.inp",'w')
b = np.zeros(tot)
c = 0

for i in range (ran):
	a = random.randrange(-5,0,4)
	for j in range(dur):
		b[c] = a
		c = c + 1

np.savetxt(f1,b)
