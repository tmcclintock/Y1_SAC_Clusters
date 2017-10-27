"""
Find the mean richnesses and redshifts of the subsamples.
"""
import numpy as np

fname = "clusters_z%0d_l%0d.txt"
Nl,Nz = 7,3
lams = np.zeros((Nz,Nl))
zs = np.zeros_like(lams)

for i in range(0,Nz):
    for j in range(0,Nl):
        data = np.loadtxt(fname%(i,j))
        lams[i,j] = np.mean(data[:,1])
        zs[i,j] = np.mean(data[:,0])

np.savetxt("mean_richnesses.txt",lams)
np.savetxt("mean_redshifts.txt",zs)
