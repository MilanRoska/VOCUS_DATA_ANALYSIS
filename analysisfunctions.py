# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:42:46 2024

@author: m.roska
"""
#%%packages
import pandas as pd
import h5py
import numpy as np
from datetime import datetime, timedelta


#%%Support Fnctions
def linearFunc(x,intercept,slope):
    y = intercept + slope*x
    return y

def convert_ldap_timestamp(ldap_timestamp):
    return pd.to_datetime(ldap_timestamp / 1e7 - 11644473600, unit='s')

#%%Main Function for loading hdf file

def load_hdf_file_VOCUS_processed(input_file_path):
    #open file
    file = h5py.File(input_file_path, 'r')
    #open folder
    peak_data_group = file['PeakData']
    timing_data_group = file['TimingData']
    #open datasets
    peak_data = peak_data_group['PeakData'][:]
    peak_table = peak_data_group['PeakTable'][:]
    buf_times = timing_data_group['BufTimes'][:]
    #open attributes
    time_stamp_start = convert_ldap_timestamp(timing_data_group.attrs['AcquisitionTimeZero'])
    #reshape data
    # Get the shape of the input array
    peak_data_original_shape = peak_data.shape
    peak_data_reshaped = peak_data.reshape((-1, peak_data_original_shape[3]))
    #reshape buftimes and convert to UTC
    buf_times_reshaped_elapsed = buf_times.reshape(-1)
    buf_times_reshaped_UTC = [time_stamp_start + timedelta(seconds=seconds) for seconds in buf_times_reshaped_elapsed]
   
    peak_table_labels = peak_table['label'].tolist()
    peak_table_labels = [byte.decode('utf-8') for byte in peak_table_labels]
    # Create a new list with "compound_" prefix added to each entry
    peak_table_labels = ['compound_' + label for label in peak_table_labels]
    #generate data frame
    data = pd.DataFrame(data = peak_data_reshaped, index = buf_times_reshaped_UTC, columns = peak_table_labels)
    #add meta data
    #add elapsed time
    data['meta_elapsed_time'] = buf_times_reshaped_elapsed
    return data



