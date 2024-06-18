# -*- coding: utf-8 -*-
"""
Created on Wed May  3 14:29:59 2023

@author: gdub5
"""

import matplotlib.pyplot as plt
import numpy as np
from CEA_LSP_polynomials import CEA_LSP_polynomials

R = 8314.462  # Universal gas constant [J/kmol-K]
T_ref = 298.15  # K
p_ref = 100000  # 1 bar = 1e5 Pa


def entropy_bar_i(compound, T, p_i):
    s_bar_0_i = CEA_LSP_polynomials(compound, T)[2]*R
    return s_bar_0_i - R*np.log(p_i/p_ref, where=p_i > 0)


def gibbs_bar_i(compound, T, p_i):
    h_bar_i = CEA_LSP_polynomials(compound, T)[1]*R*T  # Enthalpy
    s_bar_i = entropy_bar_i(compound, T, p_i)  # Entropy
    return h_bar_i - T*s_bar_i


# 'compound' values: 'H2', 'H', 'Hp', 'e', 'Ar', 'Arp'

# test = gibbs_bar_i('H2', 8000, 10e5)
# print(test/1e6)  # Gibbs in MJ

# Draw plot
x = np.linspace(0, 1, 10000000)
T = 4000
P = np.array([0.1e5, 1e5, 10e5])  # pressure in pa

for p in P:
    n_H = 2*x
    n_H2 = 1-x
    n_tot = n_H + n_H2
    y_H = n_H/n_tot
    y_H2 = n_H2/n_tot
    p_H = y_H*p
    p_H2 = y_H2*p
    G = (n_H*gibbs_bar_i("H", T, p_H) + n_H2*gibbs_bar_i("H2", T, p_H2))/1e6
    plt.plot(x, G)  # G in MJ/kmol
    xi = x[np.argmin(G)]
    print(xi)
    print(f'Mole Fraction of H at {p/1e5} bar: {2*xi/(1+xi)}')
    print(f'Mole Fraction of H2 at {p/1e5} bar: {(1-xi)/(1+xi)}\n')

plt.xlim(0, 1)
plt.legend(P/1e5)
plt.xlabel('x (degree of dissociation)')
plt.ylabel('G (MJ)')
plt.show()
