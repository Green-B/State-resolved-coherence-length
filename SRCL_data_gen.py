import numpy as np
import pickle as pkl
import sys
import os
from time import time as timenow

(L, W, n_config) = (int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])) # L = linear size, W = maximum disorder
N = L**3 # Number of atoms

# Condense cubic lattice indices to a single index with periodicity
def Hind(nx, ny, nz):
    return (L**2)*(nz%L) + L*(ny%L) + (nx%L)

# Invert single-indexing of Hind and return tuple with cubic lattice indices
def unHind(one_index):
    nx = one_index%L
    ny = int(((one_index-nx)%(L**2))/L)
    nz = int(( one_index - L*(ny%L) - (nx%L) )/(L**2))
    return (nx, ny, nz)

# Generate pristine Hamiltonian
H0 = np.zeros((N, N))
for nx in np.arange(L):
    for ny in np.arange(L):
        for nz in np.arange(L):
            H0[ Hind(nx, ny, nz) , Hind(nx+1, ny, nz) ] = 1
            H0[ Hind(nx+1, ny, nz) , Hind(nx, ny, nz) ] = 1
            H0[ Hind(nx, ny, nz) , Hind(nx, ny+1, nz) ] = 1
            H0[ Hind(nx, ny+1, nz) , Hind(nx, ny, nz) ] = 1
            H0[ Hind(nx, ny, nz) , Hind(nx, ny, nz+1) ] = 1
            H0[ Hind(nx, ny, nz+1) , Hind(nx, ny, nz) ] = 1

# Generate periodic distance metric - the distance bewteen sites in terms of cubic lattice indices with periodicity
def periodic_dist(dist_arr):
    return np.minimum( np.abs(dist_arr) , L-np.abs(dist_arr) )

R_dist = np.zeros((N, N))
for nx in np.arange(L):
    for ny in np.arange(L):
        for nz in np.arange(L):
            for nxp in np.arange(L):
                for nyp in np.arange(L):
                    for nzp in np.arange(L):
                        R_dist[Hind(nx,ny,nz), Hind(nxp,nyp,nzp)] = np.sqrt( periodic_dist(nx-nxp)**2 + periodic_dist(ny-nyp)**2 + periodic_dist(nz-nzp)**2 )

"""
# Attempt to vectorize; works but uses cubic lattice indices rathar than single indices
(nx_arr, ny_arr, nz_arr) = (np.arange(L).reshape(L,1,1,1,1,1), np.arange(L).reshape(1,L,1,1,1,1), np.arange(L).reshape(1,1,L,1,1,1))
(nxp_arr, nyp_arr, nzp_arr) = (np.arange(L).reshape(1,1,1,L,1,1), np.arange(L).reshape(1,1,1,1,L,1), np.arange(L).reshape(1,1,1,1,1,L))
def periodic_dist(dist_arr):
    return np.minimum( np.abs(dist_arr) , L-np.abs(dist_arr) )

R_dist = np.sqrt( (periodic_dist(nx_arr-nxp_arr))**2 + (periodic_dist(ny_arr-nyp_arr))**2 + (periodic_dist(nz_arr-nzp_arr))**2 )
"""

os.chdir(os.getcwd() + r"\SRCL_data")
filename = "SRCL_data_L{}_W{}_c{}_{}.pkl".format(L, W, n_config, timenow())

# All energies and SRCLs across all configurations
all_E = np.zeros(N*n_config)
all_SRCL = np.zeros(N*n_config)
for nc in np.arange(n_config):
    print("Now on configuration {} of {}".format(nc+1, n_config))
    # Add disorder and calculate eigenpairs
    H = H0 + np.diag(np.random.uniform(low=-W, high=W, size=N))
    E, eigvecs = np.linalg.eigh(H)
    # Calculate SRCL and save data
    SRCL = np.zeros(N)
    for n in np.arange(N):
        SRCL[n] = np.einsum("i,ij,j", (np.abs(eigvecs[:,n])**2), R_dist, (np.abs(eigvecs[:,n])**2))
    all_E[nc*N : (nc+1)*N] = E
    all_SRCL[nc*N : (nc+1)*N] = SRCL

with open(filename, "wb") as file_out:
    pkl.dump([all_E, all_SRCL], file_out)
