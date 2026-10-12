import numpy as np
from scipy.linalg import hessenberg
import sys

#computes Hilbert product between to vectors in the Hilbert space
def hil_prod(a,b):
    return np.dot(np.conj(a),b)

#given an hamiltonian H and an initial state psi_0, computes the lanczos coefficients and returns the Krylov basis and Hamiltonian in it
def lanczos(H, psi_0, threshold):
    #check compatibility of dimensions
    if(H.ndim != 2):
        raise ValueError("H must be a matrix")
    
    if(H.shape[0]!=H.shape[1]):
        raise ValueError("H must be a square matrix")
        
    if(H.shape[0]!=psi_0.size):
        raise ValueError("H and psi_0 must be compatible")

    #Check hermiticity of H 
    diff = H-np.conj(H.T)
    if(diff.any()!=0):
       raise ValueError("H must be hermitian")
    
    #normalize psi_0
    norm = np.sqrt(hil_prod(psi_0,psi_0))
    psi_0 /= norm

    dimHilbert = psi_0.size
    dimKrylov = int(1)

    #Rotate the base using Householder reflection such that psi_0 is vector (1,0,0,0,...)
    one_zeros = np.zeros_like(psi_0)
    one_zeros[0] = 1.+0.j

    diff = psi_0 - one_zeros
    diff_norm = np.sqrt(hil_prod(diff,diff))
    
    if(diff_norm != 0):
        #Householder matrix
        diff /= diff_norm
        House = np.identity(dimHilbert) -2.*np.outer(diff, diff)
        H = House@H@House

    #Now that H is in proper form, use Hessenberg form to compute Hamiltonian in Krylov base

    H_Hess = np.zeros((dimHilbert,dimHilbert))
    H_Hess = hessenberg(H)

    #Lanczos coefficients
    a = np.diag(H_Hess, k=0).copy()
    b = abs(np.diag(H_Hess, k=1)).copy() #take absolute value because the algorithm sometimes returns negative values
    #print("as=",a)
    #print("bs=",b)

    for bi in b:
        if bi < threshold:
            break
        dimKrylov += 1

    a.resize(dimKrylov)
    b.resize(dimKrylov)

    temp = b[:-1]

    H_Kryl = np.diag(a.real, k=0) + np.diag(temp, k=1) + np.diag(temp, k=-1)

    return H_Kryl

#given a Hamiltonian (already assumed to be expressed in Krylov basis) returns the time-averaged complexity at infinite time
def compl_inf_time(H_Kryl):
    #check that H_Kryl is symmetric
    diff = H_Kryl-H_Kryl.T
    if diff.any()!=0:
        raise ValueError("H must be symmetric and tridiagonal")

    #diagonalization
    E, psi = np.linalg.eigh(H_Kryl)
    dim = E.size

    kappa = np.zeros(dim, dtype=complex)

    for i in range(dim):
        for j in range(dim):
            c_ij = psi[i,j]
            c_0j = psi[0,j]
            kappa[i] += (c_ij*c_0j)**2


    weights = np.linspace(0,dim-1,dim)

    complexity = np.sum(weights*kappa)
    return complexity

def get_complexity(H, psi_0, threshold):
    try:
        new_H = lanczos(H.astype(np.complex256), psi_0.astype(np.complex256), threshold)
    except ValueError as e:
        print(e)
    
    try:
        compl = compl_inf_time(new_H.astype(np.complex128))
    except ValueError as e:
        print(e)

    return new_H, compl