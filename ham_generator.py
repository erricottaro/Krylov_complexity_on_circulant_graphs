import numpy as np
import networkx as nx
import scipy as sp
import krylov as kryl


#returns adjacency matrix of a given graph
def adjacency_matrix(g):
    #it is return as a scipy sparse array (whatever that is)
    H = nx.adjacency_matrix(G=g)
    #convert to numpy array
    H=H.toarray()
    return H

#returns a circulant cycle of dimension dim and with phase phi
def circulant_cycle(dim, phi):
    alpha = np.exp(-phi*1.j)
    H = np.zeros((dim,dim),dtype=complex)
    for i in range(dim):
        H[i,(i+1)%dim] = alpha

    H = np.conj(H.T) + H

    return H

#returns a variable distributed according to rho(x)=2x in [0,1)
def distr_linear(shape = None):
    y = np.random.random(size=shape)
    return np.sqrt(y)

#returns random hermitian matrix with given dimension and matrix elements uniformly distributed in the unit disk (I'M SAMPLING WRONG)
def uniform_herm(dim):
    H = np.zeros((dim,dim), dtype=complex)
    for i in range(dim):
        for j in range(i+1,dim):
            rho = distr_linear()
            phase = np.random.uniform(0.,2.*np.pi)
            H[i,j] = rho*np.exp(phase*1.j)
    
    for i in range(dim):
        rand = np.random.uniform(-1.,1.)
        H[i,i] = rand
    H = H+np.conj(H.T)
    return H

#1D Anderson Hamiltonian with given dimension and disorder parameter W
def Anderson_H(dim, W):
    g = nx.path_graph(dim)
    H = nx.to_numpy_array(G=g)
    diag = np.random.uniform(low=-W, high=W, size=dim)

    H += np.diag(diag, k=0)

    return H

#Functions to return a circulant hermitian matrix 

#shift matrix (unitary representation of Z_n group)
def shift_matrix(dim, num):
    P = np.zeros((dim,dim))
    for i in range(dim):
        P[(i+num)%dim, i] = 1
    
    return P

#circulant hermitian matrix uniquely defined by a vector of complex numbers c with length floor(dim/2)+1
def circulant_ham(c, dim):

    #compatibility check between vector c and dim
    length = len(c)
    diff = int(dim/2) + 1 - int(length)
    
    if(diff!=0):
        raise ValueError("Vector of parameters incompatible with dimension of matrix")
    
    #first parameter must be real
    if(c[0].imag!=0):
        raise ValueError("c_0 must be real")
    
    #declare Hamiltonian variable
    H = np.zeros((dim,dim),dtype=complex)

    #introduce all the terms that for sure are not redundant
    for i in range(length-2):
            H += c[i+1]*shift_matrix(dim, i+1)

    H = H + np.conj(H.T)
    #Case even
    if(dim%2==0):
        if(c[-1].imag!=0):
            raise ValueError("c_(N/2) must be real")

        H += c[-1]*shift_matrix(dim, int(dim/2))
    #case odd
    else:
        temp = c[-1]*shift_matrix(dim, length-1) + np.conj(c[-1]*shift_matrix(dim, length))
        H += temp

    #add the identity

    H += c[0]*shift_matrix(dim, 0)

    return H

#router as defined in AVS Quantum Sci. 5, 025001 (2023)
def router(gamma=0, theta=0):
    dim = 6
    H = np.zeros((dim,dim), dtype=complex)

    sub_H = circulant_ham((gamma, np.exp(theta/3.*1.j)), dim=3)

    H[1:4,1:4] = sub_H.copy()

    H[0,1] = 1
    H[1,0] = 1

    H[2,4] = 1
    H[4,2] = 1

    H[3,5] = 1
    H[5,3] = 1
    
    return H
