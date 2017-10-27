"""
Compare errorbars for different scatters.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=16)
zlenses = np.loadtxt("/Users/tmcclintock/Data/DATA_FILES/y1_data_files/Y1_meanz.txt")
h = 0.7
Nbins = 15

def get_Rbins(zi, lj): #Mpc physical
    binmin = 0.0323
    binmax = 30.0
    Redges = np.logspace(np.log(binmin), np.log(binmax), num=Nbins+1, base=np.e)
    return (Redges[:-1]+Redges[1:])/2.

def get_errorbar(inpath):
    C = np.loadtxt(inpath)
    return np.sqrt(np.diag(C))

if __name__ == "__main__":
    pss = [0, 10, 25, 45, 60] #percent scatters
    zi, lj = 0, 6
    inpath = "tom_covariance_MLps%d_z%d_l%d.txt"
    labels = [r"$%d$",r" Scatter"]
    Rbins = get_Rbins(zi, lj)
    for ps in pss:
        err = get_errorbar(inpath%(ps, zi, lj))
        label = labels[0]%ps + r"\%" + labels[1]
        if ps != 25: plt.loglog(Rbins, err, ls='-', label=label)
        else: plt.loglog(Rbins, err, ls='-', c='k', lw=2, label=label)

    types = ["MC", "ML", "MIS"]
    labels = [r"$M-c$ only", r"$M-\lambda$ 25\% only", "Miscentering only"]
    inpath = "cov_%sonly_z0_l6.txt"
    for i, typ in zip(range(len(types)), types):
        err = get_errorbar(inpath%typ)
        plt.loglog(Rbins, err, ls=':', label=labels[i])
    plt.ylabel(r"$\delta\Delta\Sigma\ [{\rm M_\odot/pc^2}]$")
    plt.xlabel(r"$R\ [{\rm Mpc}]$")
    plt.subplots_adjust(bottom=0.15, left=0.15)
    plt.gca().grid(alpha=0.5)

    base = "/Users/tmcclintock/Data/"
    shapepath = base+"DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1clust_l%d_z%d_shapecov.dat"%(lj, zi)
    cov_shape = np.loadtxt(shapepath)
    ersh = np.sqrt(np.diag(cov_shape))
    plt.plot(Rbins, ersh, c='pink', label=r"Shape")

    
    jkpath = base+"/DATA_FILES/y1_data_files/blinded_tamas_files/full-mcal-raw_y1subtr_l6_z0_dst_cov.dat"
    cov_jk = np.loadtxt(jkpath)
    erjk = np.sqrt(np.diag(cov_jk))
    plt.scatter(Rbins, erjk, c='k', label=r"JK est.")

    plt.legend(loc = 0, fontsize=12, frameon=False)
    plt.show()
