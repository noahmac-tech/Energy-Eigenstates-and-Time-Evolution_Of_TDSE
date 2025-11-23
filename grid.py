import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

L = 100
N = 500

xi = np.linspace(-L/2,L/2,N)
dx = xi[1] - xi[0]