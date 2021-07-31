#============================
# Importing python modules  #
#============================
import os, sys
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point

#============================
# User defined functions    #
#============================
# Pressure of non-interacting gas
def fPressure(p,num,T,key):

    if key == "hadrons":

        m = mass[num] # hadron mass, obtained from hadronList.txt
        g_i = g[num] # hadron degeneracy factor, obtained from hadronList.txt
        pm_1 = pm1[num] # +1: fermi-dirac distr., -1: bose-einstein distr.; obtained from hadronList.txt

    if key == "partons":

        m = mass_Partons[num] # parton mass, obtained from "partonList" defined within this code
        g_i = g_Partons[num] # parton degeneracy factor, obtained from "partonList" defined within this code
        pm_1 = pm1_Partons[num] # +1: fermi-dirac distr., -1: bose-einstein distr.; obtained from "partonList" defined within this code

    E = np.sqrt(m**2+p**2) # particle energy
    normFactor = (g_i/(6*(np.pi**2))) # normalization factor
    
    return normFactor*(1/E)*(p**4)*np.exp(-E/T)*1/(pm_1*(pm_1+np.exp(-E/T))) # Function to be integrated

#============================
# Main function             #
#============================
# Declaration of variables
# ---------------------------
Tmin = 150 # min. temperature (in MeV)
Tmax = 190 # max. temperature (in MeV)
B = 236**4 # bag constant in (MeV^4) according to hadron spectroscopy results
lowerLimit = 0 # lower limit for integral: 0
upperLimit = np.inf # upper limit for integral: +infinity
# -----------------------------
temp = []
pressureHadronGas = []
pressureQGP = []
pressureQGPMass = []
# -----------------------------

# Path to the directory containing the hadronList.txt file
path_hadrons = os.path.join(sys.path[0], "hadronList.txt")

# Reads the txt input file
hadronList = pd.read_table(path_hadrons, sep="\s+")

# Lists corresponding to columns with the particle name, g-Factor, +-1 referring to fermions/bosons and particle mass
hadrons = hadronList["col1"].values.tolist() # name of particle
g = hadronList["col3"].values.tolist() # degeneracy factor
pm1 = hadronList["col4"].values.tolist() # +1 or -1 (fermion or boson)
mGeV = hadronList["col5"].values.tolist() # particle mass (in GeV)
mass = [element * 1000 for element in mGeV] # uncomment for result in MeV

# Table for massless quarks and gluons is defined below
partonList = [["massless_quarks", "gluons", "s_quark"],[21, 16, 21/2],[1,-1,1],[0,0,95]] # for more details about g-Factor, see project
# Lists corresponding to columns with the particle name, g-Factor, +-1 referring to fermions/bosons and particle mass
partons = partonList[0] # name of particle
g_Partons = partonList[1] # degeneracy factor
pm1_Partons = partonList[2] # +1 or -1 (fermion or boson)
mass_Partons = partonList[3] # particle mass (already in MeV)

for T in np.arange(Tmin, Tmax+1, 1):

    p_HG = 0.0 # total pressure of non-interacting hadron gas
    p_QG = 0.0 # total pressure of non-interacting quarks and gluons, massless u,d quarks
    p_QG_mass = 0.0 # total pressure of non-interacting quarks and gluons, massless u, d quarks and massive s quark

    for i in range(0,len(hadrons)):

        key = "hadrons"
        # Individual and total pressure: 
        p_i, err = quad(fPressure, lowerLimit, upperLimit, args=(i, T, key)) # integral of user defined function, individual pressure
        p_HG += p_i # sum of all individual pressures, obtains total pressure of non-interacting hadron gas

    for j in range(0,len(partons)):

        key = "partons"
        p_j, err = quad(fPressure, lowerLimit, upperLimit, args=(j, T, key))  # integral of user defined function, individual pressure

        if j < 2:
            p_QG += p_j  # sum of individual pressures (massless u,d quarks; gluons), obtains total pressure of non-interacting quarks and gluons (approx) inside bag
        if j < 3:
            p_QG_mass += p_j # sum of individual pressures (massless u,d quarks, massive s quark; gluons), obtains total pressure of non-interacting quarks and gluons (realistic case) inside bag

    p_QG = p_QG-B # bag pressure subtracted from total pressure
    p_QG_mass = p_QG_mass - B # bag pressure subtracted from total pressure

    temp.append(T)
    pressureHadronGas.append(p_HG)
    pressureQGP.append(p_QG)
    pressureQGPMass.append(p_QG_mass)

#============================
# Plotting                  #
#============================

# Data plotting
plt.plot(temp,pressureHadronGas, "k", label="Hadron Gas")
plt.plot(temp,pressureQGP, "b", label="QGP; u,d quarks")
plt.plot(temp,pressureQGPMass, "g", label="QGP; u,d,s quarks")

# Axis labels and title
plt.title("Transition temperature: Hadron Gas to QGP")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (MeV$^{4}$)")

# Calculates and plots the intersection of the curves
curve_1 = LineString(np.column_stack((temp,pressureHadronGas)))
curve_2 = LineString(np.column_stack((temp,pressureQGP)))
curve_3 = LineString(np.column_stack((temp,pressureQGPMass)))

T_tr_massless = curve_1.intersection(curve_2)
T_tr_massive = curve_1.intersection(curve_3)

plt.plot(*T_tr_massless.xy, "rv")
plt.plot(*T_tr_massive.xy, "r*")
plt.text(T_tr_massless.x+0.5, T_tr_massless.y-2e8, "$T_{tr}= $" + "{:.1f}".format(T_tr_massless.x))
plt.text(T_tr_massive.x+0.5, T_tr_massive.y-2e8, "$T_{tr}= $" + "{:.1f}".format(T_tr_massive.x))

plt.legend()
plt.savefig("pressureVsTemp.png", dpi=300)
plt.show()
