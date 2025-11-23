import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# Potential Well parameters

def q_parameter(p):
    return 1 / np.sqrt(p)

# Potential Energy Function

def potential_energy(x,q):
    return - 1 / (np.cosh(q * x) ** 2)