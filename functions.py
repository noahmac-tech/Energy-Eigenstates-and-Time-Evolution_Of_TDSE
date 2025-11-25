import numpy as np
from numpy.linalg import eigh
from scipy import constants
 
def q_parameter(p):
    return 1 / np.sqrt(p)

def potential_energy(x,q):
    return - 1 / (np.cosh(q * x) ** 2)

def number_of_localized_states(p):
    s = 0.5 * (np.sqrt(1 + 4 * p) - 1)
    return int(s)

def localized_analytical_energies(p,n):
    s = number_of_localized_states(p)
    energies_analytical = []
    for n in range(int(s)):
        E_n = -(1/(4*p)) * (-(1 + 2*n) + np.sqrt(1 + 4*p))**2
        if E_n < 0: 
            energies_analytical.append(E_n)
    return energies_analytical

def stationary_hamiltonian(N,L,p,dx,xi):

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

def Normalising_wavefunction(bound_wavefunctions, n, dx):

    u = bound_wavefunctions[:, n]
    u_norm = u / np.sqrt(np.trapezoid(np.abs(u)**2, dx=dx))
    return u_norm

def Selecting_energy(bound_energies, n):

    epsilon_n = bound_energies[n]
    return epsilon_n

def Crank_Nicholson_Matrices(N,dt,H):

    I = np.identity(N)
    A = I + 0.5j *dt * H
    B = I - 0.5j *dt * H

    A_inv = np.linalg.inv(A)
    return A_inv, B

def Crank_Nicholson_Step(psi_t, A_inv, B):

    psi_t_dt = A_inv @ (B @ psi_t)
    return psi_t_dt

def Modulation_Frequency(epsilon_0, epsilon_2):

    omega = np.abs(epsilon_2 - epsilon_0)
    return omega

def Time_Evolving_Potential(xi, p, t, eta, omega):
    
    V_t = - (1 + eta * np.sin(omega * t)) / (np.cosh(q_parameter(p) * xi) ** 2)
    return V_t

def Time_Evolving_Hamiltonian(N,L,p,dx,xi,t, eta):

    V_t = Time_Evolving_Potential(xi, p, t, eta)
    main_diag = 2 / dx**2 + V_t
    off_diag = -1 / dx**2 * np.ones(N-1)

    H_t = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)
    return H_t