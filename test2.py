import numpy as np
import matplotlib.pyplot as plt

# Constants
c = 3e8  # Speed of light (m/s)
me = 9.10938356e-31  # Electron mass (kg)
re = 2.8179403227e-15  # Classical electron radius (m)
Na = 6.02214076e23  # Avogadro's number (mol^-1)
e_charge = 1.60217662e-19  # Elementary charge (C)

# Material properties of water
Z = 10  # Atomic number of water (approx. Z = 10 for H2O)
A = 18e-3  # Molar mass of water (kg/mol)
rho = 1.0e3  # Density of water (kg/m^3)
I = 75e-6 * e_charge  # Mean excitation energy of water (J)

# Initial electron energy
E_initial = 10e6 * e_charge  # 10 MeV in joules


def stopping_power_soft(E, beta, gamma):
    Tmax = 2 * me * c**2 * beta**2 * gamma**2 / (1 + 2 * gamma * me / me)  # Maximum energy transfer
    return (4 * np.pi * re**2 * me * c**2 * Z * Na / (A * beta**2)) * (np.log(Tmax / I) - beta**2)


def hard_collision_loss(beta, gamma):
    Tmax = 2 * me * c**2 * beta**2 * gamma**2 / (1 + 2 * gamma * me / me)  # Max energy transfer
    return np.random.uniform(0, Tmax)  # Random energy transfer

    # Simulation parameters
dx = 1e-7  # Step size in meters
x_max = 10  # Maximum depth in meters
num_steps = int(x_max / dx)

# Initialize arrays for depth, energy, and energy loss
depth = np.zeros(num_steps)
energy = np.zeros(num_steps)
energy_loss = np.zeros(num_steps)

# Initial conditions
energy[0] = E_initial
current_depth = 0

for i in range(1, num_steps):
    if energy[i-1] <= 0:
        break

    # Calculate velocity and relativistic factors
    gamma = energy[i-1] / (me * c**2) + 1
    beta = np.sqrt(1 - 1 / gamma**2)

    # Energy loss due to soft collisions
    dE_soft = stopping_power_soft(energy[i-1], beta, gamma) * dx

    # Energy loss due to hard collisions
    dE_hard = hard_collision_loss(beta, gamma)

    # Total energy loss for this step
    dE_total = dE_soft + dE_hard

    # Update energy and depth
    energy[i] = energy[i-1] - dE_total
    depth[i] = depth[i-1] + dx
    energy_loss[i] = dE_total



# Truncate arrays to the stopping point
valid_indices = energy > 0
depth = depth[valid_indices]
energy = energy[valid_indices]
energy_loss = energy_loss[valid_indices]

# Plot Energy vs Depth
plt.figure(figsize=(10, 6))
plt.plot(depth * 1e2, energy / e_charge, label="Energy (MeV)")
plt.xlabel("Depth (cm)")
plt.ylabel("Energy (MeV)")
plt.title("Energy Loss of 10 MeV Electron in Water")
plt.grid()
plt.legend()
plt.show()

# Plot Energy Loss Per Step
plt.figure(figsize=(10, 6))
plt.plot(depth * 1e2, energy_loss / e_charge, label="Energy Loss per Step (MeV)")
plt.xlabel("Depth (cm)")
plt.ylabel("Energy Loss per Step (MeV)")
plt.title("Energy Loss per Step in Water")
plt.grid()
plt.legend()
plt.show()
