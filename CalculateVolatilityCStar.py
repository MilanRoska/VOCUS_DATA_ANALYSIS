# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:39:24 2024

@author: m.roska

based on code by Alexandra Tsimpidi
"""

# %% packages
import re
import numpy as np
from collections import defaultdict

# %% suport functions

# count atom number
def parse_chemical_formula(formula):
    # Regex pattern to capture elements and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    # Find all matches
    matches = re.findall(pattern, formula)
    # Dictionary to store the counts of each element
    element_counts = defaultdict(int)
    for (element, count) in matches:
        if count == '':
            count = 1
        else:
            count = int(count)
        element_counts[element] += count
    return element_counts

# apply count number of atoms
def atom_num(formula, atom):
    element_counts = parse_chemical_formula(formula)
    return element_counts[atom]

# %% main function

def calc_vol(spc_eqn):
    # checks if carbon and hydrogen are in formula
    if all(element in spc_eqn for element in ['C', 'H']):
        # checks if Br Cl I and F are not in formula
        if all(element not in spc_eqn for element in ['Br', 'Cl', 'I', 'F']):
           # calcualte number of atoms
           nC = atom_num(spc_eqn,'C')
           nH = atom_num(spc_eqn,'H')
           nO = atom_num(spc_eqn,'O')
           nN = atom_num(spc_eqn,'N')
           nS = atom_num(spc_eqn,'S')
           # Apply equation of Li et al. 2016 (Molecular corridors and parameterizations of volatility in the chemical evolution of organic aerosols, ACP, 2016)
           # log10 C0 - 
           n0C = 23.80
           bC  = 0.4861
           bO  = 0.0
           bCO = 0.0
           bN  = 0.0
           bS  = 0.0
           # choose set of variables based on atom types contained
           if 'O' in spc_eqn :
              n0C = 22.66
              bC  = 0.4481
              bO  = 1.656
              bCO = -0.7790
              bN  = 0.0
              bS  = 0.0
           if 'N' in spc_eqn :
              n0C = 24.59
              bC  = 0.4066
              bO  = 0.0
              bCO = 0.0
              bN  = 0.9619
              bS  = 0.0
           if 'O' in spc_eqn and 'N' in spc_eqn :
              n0C = 24.13
              bC  = 0.3667
              bO  = 0.7732
              bCO = -0.07790
              bN  = 1.114
              bS  = 0.0
           if 'O' in spc_eqn and 'S' in spc_eqn:
              n0C = 24.06
              bC  = 0.3637
              bO  = 1.327
              bCO = -0.3988
              bN  = 0.0
              bS  = 0.7579
           if 'O' in spc_eqn and 'N' in spc_eqn and 'S' in spc_eqn:
              n0C = 28.50
              bC  = 0.3848
              bO  = 1.011
              bCO = 0.2921
              bN  = 1.053
              bS  = 1.316
           # calculate voaltility as log(C0) of formula
           vol = ((n0C-nC)*bC-nO*bO-2*(nC*nO/(nC+nO))*bCO-nN*bN-nS*bS)
        else:
            print('Br, Cl, I or F in forumla')
            vol = np.nan
    else:
        print('no C and/or H in forumla')
        vol = np.nan
    return vol

# %% run example
spc_eqn = 'C10H30O5'
vol = calc_vol(spc_eqn)