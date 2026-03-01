import numpy as np
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
    dimKrylov = int(0)

    #Lanczos coefficients
    a = np.ndarray(dimHilbert+1, dtype=complex)
    b = np.ndarray(dimHilbert+1, dtype=complex)

    #Krylov base
    Kryl = np.ndarray(shape=(dimHilbert+1, dimHilbert), dtype=complex)
    Kryl[0] = psi_0.copy()
    #print("0th Krylov element:", Kryl[0])

    #initalization
    b[0]=0
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
    H_Kryl = np.zeros(shape=(dimKrylov,dimKrylov), dtype=float)
    #set transition energies
    temp = b[1:].real
    H_Kryl = np.diag(temp,1)
    #make hermitian
    H_Kryl = H_Kryl + np.conj(H_Kryl.T)
    #set onsite energies
    np.fill_diagonal(H_Kryl, val=a.real)
    

    return Kryl, H_Kryl

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

    kappa = np.zeros(dim)

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
        base, new_H = lanczos(H, psi_0)
    except ValueError as e:
        print(e)
    
    try:
        compl = compl_inf_time(new_H)
    except ValueError as e:
        print(e)

    return new_H, compl