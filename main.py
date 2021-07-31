#============================================================================
# Importing modules
import os, sys
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pylab as plt

#===========================================================================
# Functions defined by user, compute the integral of each particle species #
#=========================================================================== 
# Hadrons
def fHG(p,i,T):
    E = np.sqrt(m[i]**2+p**2)
    normFactor = (g[i]/(6*(np.pi**2)))

    return normFactor*(1/E)*(p**4)*np.exp(-E/T)*1/(pm1[i]*(pm1[i]+np.exp(-E/T)))

# Partons
def fQG(p,j,T):
    normFactor = (g_Partons[j]/(6*(np.pi**2)))

    return normFactor*(p**3)*np.exp(-p/T)*1/(pm1_Partons[j]*(pm1_Partons[j]+np.exp(-p/T)))

# Massive s-quark
def fQG_mass(p,j,T):
    E = np.sqrt(m_Partons[j]**2+p**2)
    normFactor = (g_Partons[j]/(6*(np.pi**2)))

    return normFactor*(1/E)*(p**4)*np.exp(-E/T)*1/(pm1_Partons[j]*(pm1_Partons[j]+np.exp(-E/T)))

#============================================================================ 
# Main function                                                             #
#============================================================================
# Declaration of variables
# -----------------------------
lowerLimit = 0 # lower limit for integral: 0
upperLimit = np.inf # upper limit for integral: +infinity
Tmin = 140 # min. temperature (in MeV)
Tmax = 220 # max. temperature (in MeV)
B = 236**4 # bag constant in (MeV^4) according to hadron spectroscopy results

temp = []
pressureHadronGas = []
pressureQGP = []
pressureQGPMass = []
# -----------------------------

# Path to the directory containing the hadron list txt file
path_hadrons = os.path.join(sys.path[0], 'hadronList.txt')

# Reading the txt input file with the particle list:
hadronList = pd.read_table(path_hadrons, sep="\s+")

# Interpreting the columns containing the particle name, g-Factor, +-1 referring to fermions/bosons and particle mass
hadrons = hadronList['col1'].values.tolist() # name of particle
g = hadronList['col3'].values.tolist() # degeneracy factor
pm1 = hadronList['col4'].values.tolist() # +1 or -1 (fermion or boson)
mGeV = hadronList['col5'].values.tolist() # particle mass (in GeV)
m = [element * 1000 for element in mGeV] # uncomment for result in MeV

# Table for massless quarks and gluons is defined below
partonList = [["massless_quarks", "gluons", "s_quark"],[21, 16,21/2],[1,-1,1],[0,0,95]] #21=7/8*24

partons = partonList[0] # name of particle
g_Partons = partonList[1] # degeneracy factor
pm1_Partons = partonList[2] # +1 or -1 (fermion or boson)
m_Partons = partonList[3] # particle mass (already in MeV)

for T in np.arange(Tmin, Tmax+1, 1):

    p_HG = 0.0 # total pressure of non-interacting hadron gas
    p_QG = 0.0 # total pressure of non-interacting quarks and gluons
    p_QG_mass = 0.0 # total pressure of non-interacting quarks and gluons

    for i in range(0,len(hadrons)):

        p_i, err = quad(fHG, lowerLimit, upperLimit, args=(i, T)) # computes the integral -> pressure of each particle species
        p_HG += p_i # sum of integral results of each particle, obtains total pressure of non-interacting hadron gas

    for j in range(0,len(partons)-1):

        p_j, err = quad(fQG, lowerLimit, upperLimit, args=(j, T)) # computes the integral -> pressure of each particle species
        p_QG += p_j # sum of integral results of each particle, obtains total pressure of non-interacting hadron gas 

    p_j_mass, err_mass = quad(fQG_mass, lowerLimit, upperLimit, args=(2, T)) # computes the integral -> pressure of each particle species
    p_QG_mass += p_j_mass # sum of integral results of each particle, obtains total pressure of non-interacting hadron gas 

    p_QG = p_QG-B
    p_QG_mass = p_QG_mass + p_QG - B
    # p_QG_mass = p_QG_mass-B

    temp.append(T)
    pressureHadronGas.append(p_HG)
    pressureQGP.append(p_QG)
    pressureQGPMass.append(p_QG_mass)

plt.plot(temp,pressureHadronGas, "r", label='Hadron Gas')
plt.plot(temp,pressureQGP, "b", label='QGP, massless u,d quarks')
plt.plot(temp,pressureQGPMass, "g", label='QGP, massive s-quark')

plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (MeV$^{4}$)')

plt.title('Transition temperature of Hadron Gas to QGP')
plt.legend()
plt.show()
plt.savefig("pressureVsTemp.png")
