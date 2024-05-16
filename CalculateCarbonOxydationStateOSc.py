# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:39:24 2024

@author: m.roska
"""

# %% packages
import re
from collections import defaultdict

# %% support functions

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

# %% main functions

#calclate carbon oxidational state OSc
def calc_osc(spc_eqn):
    # calcualte number of atoms
    nC = atom_num(spc_eqn,'C')
    nH = atom_num(spc_eqn,'H')
    nO = atom_num(spc_eqn,'O')
    # calcualte Carbon Oxdational State
    OSc = (2*nO/nC) - nH/nC
    return OSc

#calculate Oxygen to Carbon ratio O:C
def calc_o_to_c(spc_eqn):
    # calcualte number of atoms
    nC = atom_num(spc_eqn,'C')
    nO = atom_num(spc_eqn,'O')
    # calcualte Oxygen to Carbon ratio
    o_to_c = nO/nC
    return o_to_c

#calculate Hydrogen to Carbon ratio H:C
def calc_h_to_c(spc_eqn):
    # calcualte number of atoms
    nC = atom_num(spc_eqn,'C')
    nH = atom_num(spc_eqn,'H')
    # calcualte Oxygen to Carbon ratio
    h_to_c = nH/nC
    return h_to_c

# %% run example
spc_eqn = 'CH4'
osc = calc_osc(spc_eqn)
o_to_c = calc_o_to_c(spc_eqn)

    