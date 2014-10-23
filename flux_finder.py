import os
import sys
import random
import numpy as np
from scipy import integrate
from scipy import constants
from scipy import fftpack
import cmath
import h5py as h5
import matplotlib.pyplot as plt

E_incd = []
E_trans_incd = []
t = []
E_refl = []
E_trans = []
E_through_incd =[]
E_through =[]

f_incd = open('/home/tap620/git/FDTD/output_data/incd/dtc_out_field_1.dat','r')
for line in f_incd:
	E_incd.append(float(line.split('\t')[3]) + 1j*float(line.split('\t')[4]))
	t.append(float(line.split('\t')[0]))
f_incd.close()
f_refl  = open('/home/tap620/git/FDTD/output_data/blocked/dtc_out_field_1.dat','r')
for line in f_refl:
	E_refl.append(float(line.split('\t')[3]) + 1j*float(line.split('\t')[4]))
f_refl.close()

f_trans = open('/home/tap620/git/FDTD/output_data/blocked/dtc_out_field_0.dat','r')
for line in f_trans:
	E_trans.append(float(line.split('\t')[3]) + 1j*float(line.split('\t')[4]))
f_trans.close()

f_trans_incd = open('/home/tap620/git/FDTD/output_data/incd/dtc_out_field_0.dat','r')
for line in f_trans_incd:
	E_trans_incd.append(float(line.split('\t')[3]) + 1j*float(line.split('\t')[4]))
f_trans_incd.close()


E_incd = np.array(E_incd)
E_refl = np.array(E_refl)

E_trans = np.array(E_trans)
E_trans_incd = np.array(E_trans_incd)


E_refl = E_refl-E_incd
#E_refl = E_refl-E_incd

E_incd_freq = fftpack.fft(E_incd)
E_refl_freq = fftpack.fft(E_refl)

E_trans_incd_freq = fftpack.fft(E_trans_incd)
E_trans_freq = fftpack.fft(E_trans)

pow_incd = E_incd_freq * E_incd_freq.conjugate() / 2.0
pow_refl = E_refl_freq * E_refl_freq.conjugate() / 2.0

pow_trans_incd = E_trans_incd_freq * E_trans_incd_freq.conjugate()/2.0
pow_trans = E_trans_freq * E_trans_freq.conjugate() / 2.0


R = pow_refl / pow_incd
T = pow_trans / pow_incd
I = pow_incd / pow_trans_incd
freq = -fftpack.fftfreq(len(R), d = 0.005)


# plt.plot(freq,E_incd_freq)

# outfile = open("RT.dat","w")
# for ii in range(len(R)):
# 	if(freq[ii] >0 and freq[ii] < 2.0):
# 		outfile.write(str(freq[ii]) + "\t\t" + str(R[ii].real) + "\t\t" + str(T[ii].real) + '\n')
# 		# outfile.write(str(freq[ii]) + "\t\t" + str(R[ii].real) + '\t\t' + str(T[ii].real) + "\t\t" + str(Th[ii].real) + '\n')
# outfile.close()


# f_E1 = h5.File('/home/tap620/meep_test/meep_test-E_1_incd.h5', 'r');
# E1 = f_E1.values()[1]
# plt.plot(E1)

f_block = open('/home/tap620/meep_test/block.dat','r')
f_incd = open('/home/tap620/meep_test/incd.dat', 'r')

T_meep = []
R_meep = []
T_incd_meep = []
R_incd_meep = []
freq_meep = []
for line in f_block:
	T_meep.append(float(line.split(',')[2]))
	R_meep.append(float(line.split(',')[3]))
	freq_meep.append(float(line.split(',')[1]))
for line in f_incd:
	T_incd_meep.append(float(line.split(',')[2]))
	R_incd_meep.append(float(line.split(',')[3]))
T_meep = np.array(T_meep)
T_incd_meep = np.array(T_incd_meep)
R_meep = np.array(R_meep)
R_incd_meep = np.array(R_incd_meep)
freq_meep = np.array(freq_meep)

# f_block_100 = open('/home/tap620/meep_test/100_block.dat','r')
# f_incd_100 = open('/home/tap620/meep_test/100_incd.dat', 'r')

# T_meep_100 = []
# R_meep_100 = []
# T_incd_meep_100 = []
# R_incd_meep_100 = []
# freq_meep_100 = []
# for line in f_block_100:
# 	T_meep_100.append(float(line.split(',')[2]))
# 	R_meep_100.append(float(line.split(',')[3]))
# 	freq_meep_100.append(float(line.split(',')[1]))
# for line in f_incd_100:
# 	T_incd_meep_100.append(float(line.split(',')[2]))
# 	R_incd_meep_100.append(float(line.split(',')[3]))
# T_meep_100 = np.array(T_meep_100)
# T_incd_meep_100 = np.array(T_incd_meep_100)
# R_meep_100 = np.array(R_meep_100)
# R_incd_meep_100 = np.array(R_incd_meep_100)
# freq_meep_100 = np.array(freq_meep_100)
plt.plot(freq_meep,R_meep/R_incd_meep, label='R_meep')
plt.plot(freq_meep, T_meep/T_incd_meep, label='T_meep')
plt.plot(freq_meep,R_meep/R_incd_meep+T_meep/T_incd_meep, label='sum_meep')

plt.plot(freq[len(freq)/2:],R[len(freq)/2:], label = 'R')
plt.plot(freq[len(freq)/2:],T[len(freq)/2:], label = 'T')
plt.plot(freq[len(freq)/2:],(R[len(freq)/2:]+T[len(freq)/2:]), label = 'sum')



# plt.plot(freq_meep,-R_meep_100/R_incd_meep_100, label='R_meep_100')
# plt.plot(freq_meep_100, T_meep_100/T_incd_meep_100, label='T_meep_100')
# plt.plot(freq_meep_100,-R_meep_100/R_incd_meep_100+T_meep_100/T_incd_meep_100, label='sum_meep_100')
# plt.plot(freq,E_through_incd_freq.real, label = 'E_through_incd_freq Re')
# plt.plot(freq,E_through_incd_freq.imag, label = 'E_through_incd_freq Im')
# plt.plot(freq,E_through_freq.real, label = 'E_through_freq Re')
# plt.plot(freq,E_through_freq.imag, label = 'E_through_freq Im')

plt.axis([0.5,2.50,0,1.2])
plt.legend()

plt.show()