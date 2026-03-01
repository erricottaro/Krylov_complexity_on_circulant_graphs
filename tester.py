import numpy as np
import krylov as kr
import ham_generator as gen
import networkx as nx
import scipy as sp
import sys
import matplotlib as plt

def get_complexity(H, psi_0):
    try:
        base, new_H = kr.lanczos(H, psi_0)
    except ValueError as e:
        print(e)
    
    try:
        compl = kr.compl_inf_time(new_H)
    except ValueError as e:
        print(e)

    return new_H, compl

if __name__ == "__main__":
    #c_0 = 0.+0.j
    #c_1 = 1.+0.j
    #c_i = 0.+1.j
    #H = np.array([[c_0,c_1, c_1],[c_1,c_0, c_1],[c_1,c_1,c_0]])
    dim = int(sys.argv[1])
    dim_fin = int(sys.argv[-1])
    #G = nx.complete_graph(dim)
    G = nx.cycle_graph(dim)
    H = gen.adjacency_matrix(G)
    print("Initial Hamiltonian:\n",H)
    psi_0 = np.zeros(dim, dtype=complex)
    psi_0[0] = 1.+0.j

    #print(psi_0)
    new_H, compl = get_complexity(H, psi_0)

    print("H_K=\n",new_H)
    print("Time average complexity:\n",compl)

    """
    indexes = np.linspace(3, dim_fin, dim_fin-2, dtype=int)
    print(indexes)
    for i in indexes:
        H = gen.adjacency_matrix(nx.complete_graph(i))
        initial = np.zeros(i, dtype=complex)
        initial[0] = 1.0+0.j
        new_H, compl = get_complexity(H, initial)
        print("Time average complexity:\n",compl)
    """

