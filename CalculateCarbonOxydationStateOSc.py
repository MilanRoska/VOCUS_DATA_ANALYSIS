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

# drop NH4+
def modify_formula(formula: str, remove_h: int, remove_n: int) -> str:
    # Step 1: Parse the formula for individual elements and their counts
    element_pattern = r'([A-Z][a-z]?)(\d*)'
    parsed_formula = re.findall(element_pattern, formula)
    
    # Step 2: Create a dictionary to store element counts
    element_counts = {}
    
    for element, count in parsed_formula:
        element_counts[element] = int(count) if count else 1
    
    # Step 3: Remove the specified number of hydrogen (H) and nitrogen (N) atoms
    if 'H' in element_counts:
        element_counts['H'] -= remove_h
        if element_counts['H'] <= 0:
            del element_counts['H']  # Remove the element if its count drops to 0 or below
    
    if 'N' in element_counts:
        element_counts['N'] -= remove_n
        if element_counts['N'] <= 0:
            del element_counts['N']  # Remove the element if its count drops to 0 or below
    
    # Step 4: Rebuild the formula without the '+' sign
    result_formula = ''
    
    for element, count in element_counts.items():
        result_formula += element
        if count > 1:
            result_formula += str(count)
    
    return result_formula


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
spc_eqn = 'C10H30O5'
osc = calc_osc(spc_eqn)
o_to_c = calc_o_to_c(spc_eqn)

    