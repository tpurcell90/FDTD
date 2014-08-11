import os
import sys
import random
import numpy as np
from scipy import integrate
from scipy import constants
from scipy import fftpack
import cmath
import matplotlib.pyplot as plt

E_incd = [];
t = [];
E_trans = [];

f_incd = open('/home/tap620/git/FDTD/output_data/incd/dtc_out_field_1.dat','r')
for line in f_incd:
	E_incd.append(float(line.split('\t')[3]))
	t.append(float(line.split('\t')[0]))
f_incd.close()
f_trans  = open('/home/tap620/git/FDTD/output_data/blocked/dtc_out_field_1.dat','r')
for line in f_trans:
	E_trans.append(float(line.split('\t')[3]))
f_trans.close()

E_incd = np.array(E_incd)
E_trans = np.array(E_trans)

#E_trans = E_trans-E_incd

E_incd_freq = fftpack.rfft(E_incd)
E_trans_freq = fftpack.rfft(E_trans)

T = E_trans_freq / E_incd_freq

pow_incd = E_incd_freq * E_incd_freq.conjugate() / 2.0
pow_trans = E_trans_freq * E_trans_freq.conjugate() / 2.0

freq = fftpack.rfftfreq(len(T), d = 0.0125)

plt.plot(T[0:100])

plt.show()