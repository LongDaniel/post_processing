#import library
import numpy as np
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import os

#Extracted individual wind turbine power
def power(i):
    t = []
    for j in range (i, int(len(x)), nturbine):
            t.append(y[j])
    return t

#-------------------------------#
#Extracted wind turbine power density
def powerdensity(i):
    powerdensity = []
    for j in range (i, int(len(x)), nturbine):
            powerdensity.append(y[j]/(Sx*Sy*D))
    return powerdensity

# create a class for turbine where store the power density and power and serial number of each turbine
class Turbine:
    def __init__(self, i = 0):
        self.i = i
        self.powerdensity = powerdensity(i)
        self.power = power(i)
    def __repr__(self):
        return 'Turbine %d' % self.i
    def __str__(self):
        return 'Turbine %d' % self.i 


#Declare working directory
os.chdir('d:\post\SWAY_9')
#Declare some variable
dt = 0.68543297937 
nturbinex = 4
nturbiney = 4
nturbine = 16
#Diameter of turbine
D = 80
#Other variable
Sx = 7
Sy = 7
#Power density
Power_density = []
#Time variable
time = []
#Mean velocity
U = 11.5258407161
#frequency of wave
omega = 0.0498561155567
#wave period
T = 2*np.pi/omega
#turbine hub height
H_hub = 70.0
#Rotational angular period
T_turb = 42.84

#Clean data
fid = open('log.uref','r')
content=fid.readlines()
# Turbine_           1 :angvel=   50.735637598099999      , TSR=   4.4059452767309510      , Uref=  0.86364504796591524
196034763     

nline = len(content)
nt = nline / nturbine
data=np.zeros((nline,4))

for i in range(int(nt)):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data[i2,0]=(i+1)*dt
    data[i2,1]=value[3] # angvel
    data[i2,2]=value[6] # TSR
    data[i2,3]=value[9] # Uref
fid.close()

np.savetxt('uref.dat', data)
fid = open('log.cthrust','r')
content=fid.readlines()
# Thrust=  0.10233803188671066      , Torque=   6.1668522808625047E-003 , C_Thrust=  0.38109671967396186      , Power=   6.3533464098887757E-002 , C_Power=  0.12840488196034763     

nline = len(content)
nt = nline / nturbine
data=np.zeros((int(nturbine*nt),6))

for i in range(int(nt)):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data[i2,0]=(i+1)*dt # time
    data[i2,1]=value[1] # thrust
    data[i2,2]=value[4] # torque
    data[i2,3]=value[7] # c_thrust
    data[i2,4]=value[10] # power
    data[i2,5]=float(value[13])*1.0 # c_power
fid.close()

np.savetxt('coeff.dat', data)
#Data Frame
data = pd.read_csv('coeff.dat', sep = '\s+', header = None)
data = pd.DataFrame(data)
x = data[0]
y = data[4] #Turbine power
#rearrange series to array multidimensional
y = np.array(y)
y = y.reshape([int(nt), nturbine])

#Normalized time variable
for i in range(0, int(len(x)), nturbine):
    #time.append(x[i]*(U/H_hub))
    time.append(x[i]/T_turb)           #Using T_turb

#calculate power density
power_density = y/(Sx*Sy*D)
power_density_average = np.mean(power_density, axis = 1)


#--------------------------------#
#Tính tay đầy thô kệch
powerdensity_1 = powerdensity(0)
powerdensity_2 = powerdensity(1)
powerdensity_3 = powerdensity(2)
powerdensity_4 = powerdensity(3)
powerdensity_5 = powerdensity(4)
powerdensity_6 = powerdensity(5)
powerdensity_7 = powerdensity(6)
powerdensity_8 = powerdensity(7)
powerdensity_9 = powerdensity(8)
powerdensity_10 = powerdensity(9)
powerdensity_11 = powerdensity(10)
powerdensity_12 = powerdensity(11)
powerdensity_13 = powerdensity(12)
powerdensity_14 = powerdensity(13)
powerdensity_15 = powerdensity(14)
powerdensity_16 = powerdensity(15)

#power density average
powerdensity_average = []
for i in range (int(len(time))):
    powerdensity_average.append((powerdensity_1[i]+ powerdensity_2[i]+powerdensity_3[i]+\
    powerdensity_4[i]+powerdensity_5[i]+powerdensity_6[i]+powerdensity_7[i]+powerdensity_8[i]+\
    powerdensity_9[i]+powerdensity_10[i]+powerdensity_11[i]+powerdensity_12[i]+powerdensity_13[i]+\
    powerdensity_14[i]+powerdensity_15[i]+powerdensity_16[i])/16)

#Normalized power density average
powerdensity_average_1 = []
for i in range(int(len(time))):
    powerdensity_average_1.append(powerdensity_average[i]/(U**3))


#--------------------------------#
#Change working directory
os.chdir('d:\post\Test_220210')

#Clean data
fid = open('log.uref','r')
content=fid.readlines()
# Turbine_           1 :angvel=   50.735637598099999      , TSR=   4.4059452767309510      , Uref=  0.86364504796591524
196034763     

nline = len(content)
nt = nline / nturbine
data=np.zeros((nline,4))

for i in range(int(nt)):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data[i2,0]=(i+1)*dt
    data[i2,1]=value[3] # angvel
    data[i2,2]=value[6] # TSR
    data[i2,3]=value[9] # Uref
fid.close()

np.savetxt('uref.dat', data)
fid = open('log.cthrust','r')
content=fid.readlines()
# Thrust=  0.10233803188671066      , Torque=   6.1668522808625047E-003 , C_Thrust=  0.38109671967396186      , Power=   6.3533464098887757E-002 , C_Power=  0.12840488196034763     

nline = len(content)
nt = nline / nturbine
data=np.zeros((int(nturbine*nt),6))

for i in range(int(nt)):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data[i2,0]=(i+1)*dt # time
    data[i2,1]=value[1] # thrust
    data[i2,2]=value[4] # torque
    data[i2,3]=value[7] # c_thrust
    data[i2,4]=value[10] # power
    data[i2,5]=float(value[13])*1.0 # c_power
fid.close()

np.savetxt('coeff.dat', data)

#plot
plt.plot(time, powerdensity_average_1, label='powerdensity_average')
plt.xlabel('(t-t0)(U/H_hub)')
plt.ylabel('Pij/U^3')
plt.title('powerdensity_average')
plt.legend()
plt.show()