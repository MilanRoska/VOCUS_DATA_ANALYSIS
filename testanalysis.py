# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:00:01 2024

@author: m.roska
"""
#%%
from analysisfunctions import load_hdf_file_VOCUS_processed

#%%
input_file_path = 'C://Users/m.roska/OneDrive - Forschungszentrum Jülich GmbH/Desktop/C3PO/Data/AllDataAmbAnd0VScan/20230330_155022_NH4_p.h5'

data = load_hdf_file_VOCUS_processed(input_file_path)
