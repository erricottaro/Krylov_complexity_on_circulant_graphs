import numpy as np

#returns a qubit in a given state on the Bloch sphere according to Haar parametrization (doesn't make sense to put here)
def qubit(theta, phi):
    c = np.cos(theta/2.)
    s = np.sin(theta/2.)
    phase = np.exp(phi*1.j)
    state = np.array([c,s*phase])
    return state

#Returns a qutrit according to Haar parametrization
def qutrit(thetas, phis):
    #rotation angles
    th1 = thetas[0]
    th2 = thetas[1]
    #complex phases
    phi1 = phis[0]
    phi2 = phis[1]

    #cosines and sines
    c1 = np.cos(th1/2.)
    s1 = np.sin(th1/2.)

    c2 = np.cos(th2/2.)
    s2 = np.sin(th2/2.)

    #phases
    phase1 = np.exp(phi1*1.j)
    phase2 = np.exp(phi2*1.j)

    #coefficients
    coeff1 = c1
    coeff2 = s1*c2*phase1
    coeff3 = s1*s2*phase2

    #state
    state = np.array([coeff1, coeff2, coeff3])

    return state
