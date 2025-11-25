import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

def Grid(L,N):

    xi = np.linspace(-L/2,L/2,N)
    dx = xi[1] - xi[0]
    return xi, dx