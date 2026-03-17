import numpy as np
from scipy.linalg import hessenberg
import sys

#computes Hilbert product between to vectors in the Hilbert space
def hil_prod(a,b):
    return np.dot(np.conj(a),b)

#given an hamiltonian H and an initial state psi_0, computes the lanczos coefficients and returns the Krylov basis and Hamiltonian in it
def lanczos(H, psi_0):
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

    threshold = 1e-6

    for bi in b:
        if bi < threshold:
            break
        dimKrylov += 1

    a.resize(dimKrylov)
    b.resize(dimKrylov)

    temp = b[:-1]

    H_Kryl = np.diag(a.real, k=0) + np.diag(temp, k=1) + np.diag(temp, k=-1)

    #Find first b_n = 0 and reshape H_kryl
    '''
    #Lanczos coefficients
    a = np.ndarray(dimHilbert+1, dtype=np.complex256)
    b = np.ndarray(dimHilbert+1, dtype=np.complex256)

    #Krylov base
    Kryl = np.ndarray(shape=(dimHilbert+1, dimHilbert), dtype=np.complex256)
    Kryl[0] = psi_0.copy()
    #print("0th Krylov element:", Kryl[0])

    #initalization
    b[0]=np.complex256(0.)
    a[0]=hil_prod(Kryl[0], H@Kryl[0])
    #print("Average energy:",a[0])
    

    #first step of algorithm
    Kryl[1] = H@Kryl[0].copy() - a[0]*Kryl[0].copy()
    b[1] = np.sqrt(hil_prod(Kryl[1],Kryl[1]))
    #print("stddv of energy:", b[1])
    i = int(1)

    #set threshold to avoid problems due to rounding errors. If b is less then such threshold, the algorithm stops (absolutely arbitrary)
    threshold = 1e-6

    #If current_b is compatible with 0, then halt the algorithm, otherwise go on
    while(b[i].real > threshold and i<dimHilbert):
        #renormalize the krylov state
        Kryl[i] /= b[i]
        #print(i,"th Krylov element:", Kryl[i])
        #compute a
        a[i] = hil_prod(Kryl[i], H@Kryl[i])

        #Lanczos step
        Kryl[i+1] = H@Kryl[i].copy() - a[i]*Kryl[i].copy() -b[i]*Kryl[i-1].copy()
        
        b[i+1] = np.sqrt(hil_prod(Kryl[i+1], Kryl[i+1]))
        #print("b[",i+1,"]:", b[i+1])

        i+=1

    #Find dimension of Krylov subspace and reshape the arrays to contain only the meaningful ones
    dimKrylov += i
    #print("dimKrylov:",dimKrylov)
    a.resize(dimKrylov)
    b.resize(dimKrylov)
    Kryl.resize((dimKrylov,dimHilbert))
    #print("as:", a)
    #print("bs:", b)
    #print(Kryl.shape)
    #print("Krylov base:", Kryl)

    #construction of Hamiltonian in Krylov base
    H_Kryl = np.zeros(shape=(dimKrylov,dimKrylov), dtype=np.complex256)
    #set transition energies
    temp = b[1:].real
    H_Kryl = np.diag(temp,1)
    #make hermitian
    H_Kryl = H_Kryl + np.conj(H_Kryl.T)
    #set onsite energies
    np.fill_diagonal(H_Kryl, val=a.real)
    '''

    return H_Kryl

#given an hamiltonian (already assumed to be expressed in Krylov basis) returns the time-averaged complexity at infinite time
def compl_inf_time(H_Kryl):
    #check that H_Kryl is symmetric
    diff = H_Kryl-H_Kryl.T
    if diff.any()!=0:
        raise ValueError("H must be symmetric and tridiagonal")

    #diagonalization
    E, psi = np.linalg.eigh(H_Kryl)
    dim = E.size
    #print("Eigenvalues:\n", E)
    #print("Eigenvectors:\n", psi)

    #print("Are eigenstates orthonormal?") #almost
    #test = (psi.T)@psi
    #print("V^T*V=",test) 
    #for i in range(psi.shape[1]):
    #    temp = psi[i,:]
    #    print(hil_prod(temp, temp))

    kappa = np.zeros(dim, dtype=complex)

    for i in range(dim):
        for j in range(dim):
            c_ij = psi[i,j]
            #print("c_",i,j,":", c_ij)
            c_0j = psi[0,j]
            #print("c_0",j,":", c_0j)
            kappa[i] += (c_ij*c_0j)**2
        #print("K_",i,":",kappa[i], "\n")


    weights = np.linspace(0,dim-1,dim)
    #print(kappa)
    #print(weights)
    #print(weights*kappa)

    complexity = np.sum(weights*kappa)
    return complexity

def get_complexity(H, psi_0):
    try:
        new_H = lanczos(H.astype(np.complex256), psi_0.astype(np.complex256))
    except ValueError as e:
        print(e)
    
    try:
        compl = compl_inf_time(new_H.astype(np.complex128))
    except ValueError as e:
        print(e)

    return new_H, compl