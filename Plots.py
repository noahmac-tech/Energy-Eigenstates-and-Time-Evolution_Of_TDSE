from Functions import *
from Grid import *
import matplotlib.pyplot as plt

# Plot 1: N = 100, L = 10, P = 30

xi, dx = Grid(L=10, N=100)
H = stationary_hamiltonian(N=100, L=10, p=30, dx=dx, xi=xi)
energies, wavefunctions = eigenvalues_and_vectors(H)
bound_energies = Bound_energies(energies)
bound_wavefunctions = Bound_wavefuncs(wavefunctions, energies)

fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(1, 1, 1)
for i, psi in enumerate(bound_wavefunctions.T):
    psi_norm = psi / np.sqrt(np.trapezoid((np.abs(psi)**2), xi))  # Normalizing the wavefunction
    ax1.plot(xi, psi_norm, label=f'n={i}, E={bound_energies[i]:.3f}')
ax1.set_xlabel('$\u03BE$')
ax1.set_ylabel('u($\u03BE$)')
ax1.legend()
plt.show()

# Plot 2: N = 100, L = 50, P = 30

xi, dx = Grid(L=50, N=100)
H = stationary_hamiltonian(N=100, L=50, p=30, dx=dx, xi=xi)
energies, wavefunctions = eigenvalues_and_vectors(H)
bound_energies = Bound_energies(energies)
bound_wavefunctions = Bound_wavefuncs(wavefunctions, energies)

fig2 = plt.figure(figsize=(10, 6))
ax2 = fig2.add_subplot(1, 1, 1)
for i, psi in enumerate(bound_wavefunctions.T):
    psi_norm = psi / np.sqrt(np.trapezoid((np.abs(psi)**2), xi))  # Normalizing the wavefunction
    ax2.plot(xi, psi_norm, label=f'n={i}, E={bound_energies[i]:.3f}')
ax2.set_xlabel('$\u03BE$')
ax2.set_ylabel('u($\u03BE$)')
ax2.legend()
plt.show()

# Plot 3: N = 100, L = 100, P = 30

xi, dx = Grid(L=100, N=100)
H = stationary_hamiltonian(N=100, L=100, p=30, dx=dx, xi=xi)
energies, wavefunctions = eigenvalues_and_vectors(H)
bound_energies = Bound_energies(energies)
bound_wavefunctions = Bound_wavefuncs(wavefunctions, energies)

fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(1, 1, 1)
for i, psi in enumerate(bound_wavefunctions.T):
    psi_norm = psi / np.sqrt(np.trapezoid((np.abs(psi)**2), xi))  # Normalizing the wavefunction
    ax1.plot(xi, psi_norm, label=f'n={i}, E={bound_energies[i]:.3f}')
ax1.set_xlabel('$\u03BE$')
ax1.set_ylabel('u($\u03BE$)')
ax1.legend()
plt.show()

# Plot 4: N = 100, L = 150, P = 30

xi, dx = Grid(L=150, N=100)
H = stationary_hamiltonian(N=100, L=150, p=30, dx=dx, xi=xi)
energies, wavefunctions = eigenvalues_and_vectors(H)
bound_energies = Bound_energies(energies)
bound_wavefunctions = Bound_wavefuncs(wavefunctions, energies)

fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(1, 1, 1)
for i, psi in enumerate(bound_wavefunctions.T):
    psi_norm = psi / np.sqrt(np.trapezoid((np.abs(psi)**2), xi))  # Normalizing the wavefunction
    ax1.plot(xi, psi_norm, label=f'n={i}, E={bound_energies[i]:.3f}')
ax1.set_xlabel('$\u03BE$')
ax1.set_ylabel('u($\u03BE$)')
ax1.legend()
plt.show()


# Question 3 

xi, dx = Grid(L=100, N=1000)
H = stationary_hamiltonian(N=1000, L=100, p=30, dx=dx, xi=xi)
energies, wavefunctions = eigenvalues_and_vectors(H)
bound_energies = Bound_energies(energies)
bound_wavefunctions = Bound_wavefuncs(wavefunctions, energies)  

u_0_norm = Normalising_wavefunction(bound_wavefunctions, 0, dx)
epsilon_0 = Selecting_energy(bound_energies, 0)

dt = 0.1
t_max = 8 * np.pi / np.abs(epsilon_0)
num_steps = int(t_max / dt)
T0 = 2 * np.pi / np.abs(epsilon_0)

c_0 = np.zeros(num_steps, dtype=complex)
time_array = np.zeros(num_steps)

for n in range(num_steps):
    
    c_0[n] = np.sum(psi * u_0_norm) * dx
    time_array[n] = n * dt
    A_inv, B = Crank_Nicholson_Matrices(N=1000, dt=dt, H=H)
    psi = Crank_Nicholson_Step(psi, A_inv, B)
