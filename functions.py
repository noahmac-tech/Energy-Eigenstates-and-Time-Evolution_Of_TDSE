import numpy as np
from numpy.linalg import eigh
from scipy import constants
from grid import xi, dx

def q_parameter(p):
    return 1 / np.sqrt(p)

def potential_energy(x,q):
    return - 1 / (np.cosh(q * x) ** 2)

def number_of_localized_states(p):
    s = 0.5 * (np.sqrt(1 + 4 / q_parameter(p)**2) - 1)
    return int(s)

def localized_analytical_energies(p,n):
    s = number_of_localized_states(p)
    energies_analytical = []
    for n in range(int(s)):
        E_n = -(1/(4*p)) * (-(1 + 2*n) + np.sqrt(1 + 4*p))**2
        if E_n < 0: 
            energies_analytical.append(E_n)
    return energies_analytical

def stationary_hamiltonian(N,L,p):

    main_diag = 2 / dx**2 + potential_energy(xi,q_parameter(p))
    off_diag = -1 / dx**2 * np.ones(N-1)

    H = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)
    return H

def eigenvalues_and_vectors(H):

    energies, wavefuncs = eigh(H)
    return energies, wavefuncs

def Bound_energies(energies):
    
    bound = np.where(energies < 0)[0]
    return energies[bound]

def Bound_wavefuncs(wavefuncs, energies):
    
    bound = np.where(energies < 0)[0]
    return wavefuncs[:, bound]

  
    
