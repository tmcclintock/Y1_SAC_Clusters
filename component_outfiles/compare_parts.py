import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=14)

zi, lj = 0, 6
mcpath ="cov_MConly_z%d_l%d.txt"
cov_mc = np.loadtxt(mcpath%(zi, lj))
mlpath ="cov_MLonly_z%d_l%d.txt"
cov_ml = np.loadtxt(mlpath%(zi, lj))
mispath ="cov_MISonly_z%d_l%d.txt"
cov_mis = np.loadtxt(mispath%(zi, lj))

ermc = np.sqrt(np.diag(cov_mc))
erml = np.sqrt(np.diag(cov_ml))
ermis = np.sqrt(np.diag(cov_mis))
erall = np.sqrt(ermc**2 + erml**2 + ermis**2)

#fullpath = "../output_files/scatter_probs/tom_covariance_z%d_l%d.txt"%(zi, lj)
#covfull = np.loadtxt(fullpath)
#erfull = np.sqrt(np.diag(covfull))

base = "/Users/tmcclintock/Data/"
shapepath = base+"DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1clust_l%d_z%d_shapecov.dat"%(lj, zi)
cov_shape = np.loadtxt(shapepath)
ersh = np.sqrt(np.diag(cov_shape))

jkpath = base+"/DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1subtr_l6_z0_dst_cov.dat"
cov_jk = np.loadtxt(jkpath)
erjk = np.sqrt(np.diag(cov_jk))


Nbins = len(ersh)
binmin = 0.0323 #Mpc phys
binmax = 30.0 #Mpc phys
Redges = np.logspace(np.log(binmin), np.log(binmax), num=Nbins+1, base=np.e)
Rbins = (Redges[:-1]+Redges[1:])/2.

x = np.arange(len(ermc))
plt.plot(Rbins, ermc, c='b', label=r"$M-c$")
plt.plot(Rbins, erml, c='r', label=r"$M-\lambda$")
plt.plot(Rbins, ermis, c='g', label=r"Mis.")
plt.plot(Rbins, ersh, c='pink', label=r"Shape")
plt.plot(Rbins, erall, c='k', ls='--', label=r"$\sum$(B,R,G)")
#plt.plot(Rbins, erfull, c='k', label=r"B,R,G all on")
plt.scatter(Rbins, erjk, c='k', label=r"JK est.")
plt.xlabel(r"$R\ [{\rm Mpc}]$")
plt.ylabel(r"$\delta\Delta\Sigma$")
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=0)
plt.title(r"$z%d\ \lambda%d$ bin"%(zi, lj))
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)
plt.gca().get_xaxis().set_ticks([])
plt.gca().get_yaxis().set_ticks([])

plt.show()
