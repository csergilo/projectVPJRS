#============================================================================
# Importing libraries
import os, sys
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt 
from math import exp, pi, sqrt


def f(x,T):
    
    return T*x**2

Tmin = 0.100 # min. temperature (GeV)
Tmax = 0.200 # max. temperature (GeV)

for T in np.arange(Tmin, Tmax, 0.010):

    print(T)
    
    integral, err = quad(f, 0, 1, args=(T)) # computes the integral -> pressure of each particle species

    print(integral)
