import numpy as np
import matplotlib.pyplot as plt
import os,sys
from os import listdir
from os.path import isfile, join

mypath = "/home/tap620/git/FDTD/fout/Ez"
onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
field_slices = []
ii = 0
for f in onlyfiles:
    print(ii)
    field_file = open('fout/Ez/'+f)
    x_val = []
    y_val = []
    field_val = []
    for line in field_file:
        x_val.append(float(line.split('\t')[0]))
        y_val.append(float(line.split('\t')[1]))
        field_val.append(float(line.split('\t')[2]))
    field_file.close()
    field = np.ndarray(shape=(max(x_val)+1,max(y_val)+1), dtype=float)
    for jj in range(len(x_val)):
        field[x_val[jj],y_val[jj]] = field_val[jj]
    field_slices.append(field)
    plt.imshow(field,extent=[-5,5,-25,25], aspect='auto')
    plt.colorbar()
    plt.savefig('fout/Ez/img/'+f[:-4] + ".png")
    plt.clf()
    ii +=1
    del field, field_val,x_val,y_val

#field_slices = np.array(field_slices)
#plt.imshow(field_slices[100,:,:],extent=[-10,10,-10,10], aspect='auto')
#plt.colorbar()
#plt.show()