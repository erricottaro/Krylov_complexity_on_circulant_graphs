import numpy as np

#returns a qubit in a given state on the Bloch sphere according to Haar parametrization
def qubit(theta, phi):
    c = np.cos(theta/2.)
    s = np.sin(theta/2.)
    phase = np.exp(phi*1.j)
    state = np.array([c,s*phase])
    return state

#returns a qubit given cartesian coordinates on a CP1
def qubit_XYZ(x):
    theta = np.ndarray(x.shape[0])
    phi = np.ndarray(x.shape[0])

    z0 = x[:,0]
    z1 = x[:,1]

    theta = 2*np.arccos(z0)
    phi = np.angle(z1)

    state = qubit(theta, phi)

    return state

#Returns a qutrit according to Haar parametrization
def qutrit(thetas, phis):
    #rotation angles
    th1 = thetas[:,0]
    th2 = thetas[:,1]
    #complex phases
    phi1 = phis[:,0]
    phi2 = phis[:,1]

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

#returns a qutrit given cartesian coordinates on a S5/S1
def qutrit_XYZ(x):
    theta1 = np.ndarray(x.shape[0])
    theta2 = np.ndarray(x.shape[0])
    phi1 = np.ndarray(x.shape[0])
    phi2 = np.ndarray(x.shape[0])

