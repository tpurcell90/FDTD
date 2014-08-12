import os
import sys
import random
import numpy as np
from scipy import integrate
from scipy import constants
from scipy import fftpack
import cmath
import matplotlib.pyplot as plt

E_incd = []
E_trans_incd = []
t = []
E_refl = []
E_trans = []

f_incd = open('/home/tap620/git/FDTD/output_data/incd/dtc_out_field_1.dat','r')
for line in f_incd:
	E_incd.append(float(line.split('\t')[3]))
	t.append(float(line.split('\t')[0]))
f_incd.close()
f_refl  = open('/home/tap620/git/FDTD/output_data/blocked/dtc_out_field_1.dat','r')
for line in f_refl:
	E_refl.append(float(line.split('\t')[3]))
f_refl.close()

f_trans = open('/home/tap620/git/FDTD/output_data/blocked/dtc_out_field_0.dat','r')
for line in f_trans:
	E_trans.append(float(line.split('\t')[3]))
f_trans.close()

f_trans_incd = open('/home/tap620/git/FDTD/output_data/incd/dtc_out_field_0.dat','r')
for line in f_trans_incd:
	E_trans_incd.append(float(line.split('\t')[3]))
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
pow_trans = E_trans_freq * E_trans_freq.conjugate() / 2.0 * 3.0

R = pow_refl / pow_incd
T = pow_trans / pow_incd
freq = fftpack.fftfreq(len(R), d = 0.01)

outfile = open("RT.dat","w")
for ii in range(len(T)):
	if(freq[ii] >0 and freq[ii] < 2.0):
		outfile.write(str(freq[ii]) + "\t\t" + str(R[ii].real) + '\t\t' + str(T[ii].real) + '\n')
outfile.close()

plt.legend()

plt.plot(freq[0:len(freq)/2],R[0:len(R)/2],freq[0:len(freq)/2],T[0:len(T)/2], freq[0:len(freq)/2.0], (R+T)[0:len(freq)/2])

plt.show()