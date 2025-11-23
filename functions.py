import numpy as np
from numpy.linalg import eigh
import matplotlib.pyplot as plt
from scipy import constants
from grid.py import xi, dx

# Potential Well parameters

def q_parameter(p):
    return 1 / np.sqrt(p)

# Potential Energy Function

def potential_energy(x,q):
    return - 1 / (np.cosh(q * x) ** 2)

def stationary_hamiltonian(N,L,p):

    main_diag = 2 / dx**2 + potential_energy(xi,q_parameter(p))
    off_diag = -1 / dx**2 * np.ones(N-1)

    H = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)
    return H

def eigenvalues_and_vectors(H):

    energies, wavefuncs = eigh(H)
    return energies, wavefuncs

