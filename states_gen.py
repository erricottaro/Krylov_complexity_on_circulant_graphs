# Krylov complexity of Chiral Discrete systems
# Copyright (C) 2026 Edoardo Ferretti
#
# This file is part of Krylov complexity of Chiral Discrete systems.
#
# Krylov complexity of Chiral Discrete systems is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# Krylov complexity of Chiral Discrete systems is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Krylov complexity of Chiral Discrete systems. If not, see <https://www.gnu.org/licenses/>.

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
    theta = np.ndarray(x.shape[0], dtype=float)
    phi = np.ndarray(x.shape[0])

    z0 = x[:,0]
    z1 = x[:,1]

    theta = 2*np.arccos(abs(z0))
    phi = np.angle(z1)

    #state = qubit(theta, phi)

    return theta, phi

#Returns a qutrit according to Haar parametrization
def qutrit(theta1, theta2, phi1, phi2):
    #cosines and sines
    c1 = np.cos(theta1/2.)
    s1 = np.sin(theta1/2.)

    c2 = np.cos(theta2/2.)
    s2 = np.sin(theta2/2.)

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

#returns a qutrit given cartesian coordinates on CP2
def qutrit_XYZ(x):
    theta1 = np.ndarray(x.shape[0])
    theta2 = np.ndarray(x.shape[0])
    phi1 = np.ndarray(x.shape[0])
    phi2 = np.ndarray(x.shape[0])

    z0=x[:,0]
    z1=x[:,1]
    z2=x[:,2]

    theta1 = 2*np.arccos(z0)
    theta2 = 2*np.arctan(abs(z2), abs(z1))

    phi1 = np.angle(z1)
    phi2 = np.angle(z2)

    #state = qutrit(theta1, theta2, phi1, phi2)

    return theta1, theta2, phi1, phi2


