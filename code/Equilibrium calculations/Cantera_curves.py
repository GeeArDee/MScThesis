# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 12:22:14 2023

@author: gdub5
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines

gas = ct.Solution("Ionized_Argon_and_Hydrogen.yaml")


# Initialize range
T_vector = np.arange(500, 20000, step=1000)
P_vector = np.array([0.1, 1, 10])*1e5

# Initialize solution vectors
Cp_vector = np.zeros((np.size(T_vector), np.size(P_vector)))
Ionization_vector = np.zeros((np.size(T_vector), np.size(P_vector)))

for i, P in enumerate(P_vector):
    # Let's do argon! (and hydrogen!)
    for j, T in enumerate(T_vector):
        gas.TPX = T, P, 'H2:1'
        gas.equilibrate('TP', solver='gibbs')
        Cp_vector[j, i] = gas.cp
        Ionization_vector[j, i] = gas.Y[2]/(gas.Y[1]+gas.Y[2])

plt.figure()
plt.plot(T_vector, Cp_vector, label=['0.1 bar', '1 bar', '10 bar'])
plt.xlabel('Temperature [K]')
plt.ylabel(r'Heat capacity $C_p$ [J/kg/K]')
labelLines()

plt.figure()
plt.plot(T_vector, Ionization_vector, label=["0.1 bar", '1 bar', '10 bar'])
plt.xlabel('Temperature [K]')
plt.ylabel(r'Degree of ionization $\alpha$ [-]')
labelLines()
