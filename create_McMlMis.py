"""
Create the part of the SAC that has M-c scatter, M-lambda scatter, and miscentering.
"""
import sys, os
import numpy as np
import helper_functions as HF
import cluster_toolkit as ct
import matplotlib.pyplot as plt

#Get concentration as a function of M and z 
concentration_spline = HF.get_concentration_spline()

#Get the cosmology dictionary
cosmo_dict = HF.get_cosmo_dict()
h  = cosmo_dict['h']
om = cosmo_dict['om']

N_realizations = 1000
N_Radii = 500

ML_percent_scatter = 25

#Power spectrum at the redshift of each clusters
#laptop:
base = "/Users/tmcclintock/Data/DATA_FILES/y1_data_files/"
#calvin
base = "/calvin1/tmcclintock/DES_DATA_FILES/DATA_FILES/y1_data_files/"
P_file_path = base+"/P_files/"
#Cluster information for each stack
cluster_file_path = base+"/cluster_files/clusters_z%d_l%d.txt"
for i in range(2, -1, -1): #z index 2, 1, 0
    for j in range(1, -1, -1): #lambda index 6 to 3, not doing 2,1,0
        #Start by getting xi_mm, which doesn't depend on mass
        k = np.loadtxt(P_file_path+"k.txt")
        Plin = np.genfromtxt(P_file_path+"./plin_z%d_l%d.txt"%(i,j))
        Pnl  = np.genfromtxt(P_file_path+"/pnl_z%d_l%d.txt"%(i,j))
        zs, lams = np.loadtxt(cluster_file_path%(i, j)).T
        zlens = np.mean(zs)
        R     = np.logspace(-2, 3, N_Radii, base=10) #go higher than BAO
        xi_mm = ct.xi.xi_mm_at_R(R, k, Pnl, exact=True)
        R_perp = np.logspace(-2, 2.4, N_Radii, base=10)

        DeltaSigma_realizations = np.zeros((N_realizations, N_Radii))
        print "Starting realizations for z%d l%d"%(i,j)
        
        for real in range(N_realizations):
            M, conc, Rmis, ismis = HF.get_cluster_parameters(lams, zs, concentration_spline, ML_scatter=ML_percent_scatter/100.)
            N_kept = len(M)
            mean_DeltaSigma = np.zeros_like(R_perp)
            biases = ct.bias.bias_at_M(M, k, Plin, om)
            
            for cl in range(N_kept): #Loop over clusters
                xi_nfw = ct.xi.xi_nfw_at_R(R, M[cl], conc[cl], om)
                bias = biases[cl]
                xi_2halo = ct.xi.xi_2halo(bias, xi_mm)
                xi_hm    = ct.xi.xi_hm(xi_nfw, xi_2halo)
                Sigma    = ct.deltasigma.Sigma_at_R(R_perp, R, xi_hm, M[cl], conc[cl], om)
                if not ismis[cl]: #isn't miscentered
                    DeltaSigma = ct.deltasigma.DeltaSigma_at_R(R_perp, R_perp, Sigma, M[cl], conc[cl], om)
                else: #is miscentered
                    Sigma_single  = ct.miscentering.Sigma_mis_single_at_R(R_perp, R_perp, Sigma, M[cl], conc[cl], om, Rmis[cl])
                    DeltaSigma = ct.miscentering.DeltaSigma_mis_at_R(R_perp, R_perp, Sigma_single)
                mean_DeltaSigma += DeltaSigma/N_kept
                
            DeltaSigma_realizations[real] = mean_DeltaSigma
            if real%100 == 0:
                print "Made realization %d"%real
        print "Made realizations for z%d l%d"%(i,j)
        #np.savetxt("output_files/stack_realizations_MLps%d_z%d_l%d.txt"%(ML_percent_scatter, i, j), DeltaSigma_realizations)
        np.savetxt("fiducial_covariances/stack_realizations/stack_realizations_z%d_l%d.txt"%(i, j), DeltaSigma_realizations)
