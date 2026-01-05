# -*- coding: utf-8 -*-
# %%Packages
from Common_Functions.StandardSetupAndCommonFunctions import convert_ldap_timestamp

import os
import io
import pandas as pd
import numpy as np
import h5py
from datetime import timedelta, datetime

# disable warning for Lines like peak_data_active_compounds_df['BufTimes'] = buf_times
pd.options.mode.chained_assignment = None  # default='warn'

# %% support functions


def load_hdf_file(input_file_full_name, compounds=None, all_compounds=True):
    print('load '+input_file_full_name)
    #load file
    current_file = h5py.File(os.path.join(input_file_full_name), 'r')
    #open groups
    peak_data_grp = current_file['PeakData']
    timing_data_grp = current_file['TimingData']
    #load data from groups
    peak_data = np.array(peak_data_grp['PeakData'])
    buf_times = np.array(timing_data_grp['BufTimes'])
    peak_table = peak_data_grp['PeakTable']
    #convert compound data to lists and decode
    formulas = peak_table['label'].tolist()[:]
    masses = peak_table['mass'].tolist()[:]
    formulas = [byte.decode('utf-8') for byte in formulas]
    
    time_stamp_0 = convert_ldap_timestamp(timing_data_grp.attrs['AcquisitionTimeZero'], unit='s')
    peak_data = peak_data.reshape((peak_data.shape[0]*peak_data.shape[1],peak_data.shape[3]))
    buf_times = buf_times.reshape((buf_times.shape[0]*buf_times.shape[1]))
    buf_times = [time_stamp_0 + timedelta(seconds=seconds) for seconds in buf_times]
    #gather Peak Data in DF
    peak_data_df = pd.DataFrame(data=peak_data, columns=formulas)
    #gather Peak Table in DF
    peak_table_df = pd.DataFrame()
    peak_table_df['formulas'] = formulas
    peak_table_df['masses'] = masses   
    # Select the columns that match the compounds array
    if all_compounds == True:
        peak_data_active_compounds_df = peak_data_df
    else:
        selected_columns = [col for col in peak_data_df.columns if col in compounds]
        # Create a new DataFrame with only the selected columns
        peak_data_active_compounds_df = peak_data_df[selected_columns]
    #check if there is an offset between buftimes and peakdata and cut of buftiems accordingly
    len_dif = len(buf_times)-len(peak_data_active_compounds_df)
    buf_times = buf_times[len_dif:]
    peak_data_active_compounds_df['BufTimes'] = buf_times
    
    # Identify and keep only the first occurrence of each unique 'BufTimes' value
    peak_data_active_compounds_df = peak_data_active_compounds_df[~peak_data_active_compounds_df['BufTimes'].duplicated(keep='first')]
    # Add the 'BufTimes' column as the index
    peak_data_active_compounds_df.set_index('BufTimes', inplace=True)
    # Convert the index to pandas datetime format
    if not isinstance(peak_data_active_compounds_df.index, pd.DatetimeIndex):
        peak_data_active_compounds_df.index = pd.to_datetime(peak_data_active_compounds_df.index)
    print(peak_data_active_compounds_df.index[0])
    print(time_stamp_0)
    return peak_data_active_compounds_df


def load_tps_from_raw_old(file_path):
    #load file
    current_file = h5py.File(file_path, 'r')
    #open groups
    tps_data_grp = current_file['TPS2']
    timing_data_grp = current_file['TimingData']
    #load data from groups
    tps_data = np.array(tps_data_grp['TwData'])
    buf_times = np.array(timing_data_grp['BufTimes'])
    tps_labels = tps_data_grp['TwInfo']
    # conert time
    time_stamp_0 = convert_ldap_timestamp(timing_data_grp.attrs['AcquisitionTimeZero'], unit='s')
    # reshape buffs
    tps_data = tps_data.reshape((tps_data.shape[0]*tps_data.shape[1],tps_data.shape[2]))
    buf_times = buf_times.reshape((buf_times.shape[0]*buf_times.shape[1]))
    # convert time
    buf_times = [time_stamp_0 + timedelta(seconds=seconds) for seconds in buf_times]
    #gather Data in DF
    tps_data_df = pd.DataFrame(data=tps_data, columns=tps_labels)
    # make buf times the index
    tps_data_df['BufTimes'] = buf_times
    tps_data_df.set_index('BufTimes', inplace=True)
    # decode labels
    tps_data_df.columns = [col.decode('latin1') if isinstance(col, bytes) else str(col) for col in tps_data_df.columns]
    # Convert the index to pandas datetime format
    if not isinstance(tps_data_df.index, pd.DatetimeIndex):
        tps_data_df.index = pd.to_datetime(tps_data_df.index)
    return tps_data_df

def load_tps_from_raw(file_path):
    import h5py
    import numpy as np
    import pandas as pd
    from datetime import timedelta

    # Load file
    current_file = h5py.File(file_path, 'r')

    # Open groups
    tps_data_grp = current_file['TPS2']
    timing_data_grp = current_file['TimingData']

    # Load data from groups
    tps_data = np.array(tps_data_grp['TwData'])
    buf_times = np.array(timing_data_grp['BufTimes'])
    tps_labels = tps_data_grp['TwInfo']

    # Convert time zero
    time_stamp_0 = convert_ldap_timestamp(timing_data_grp.attrs['AcquisitionTimeZero'], unit='s')

    # Reshape buffers
    tps_data = tps_data.reshape((tps_data.shape[0] * tps_data.shape[1], tps_data.shape[2]))
    buf_times = buf_times.reshape(-1)

    # Convert to native Python datetime list
    buf_times = [time_stamp_0 + timedelta(seconds=float(s)) for s in buf_times]

    # Gather data in DataFrame
    tps_data_df = pd.DataFrame(data=tps_data, columns=tps_labels)

    # Decode labels
    tps_data_df.columns = [col.decode('latin1') if isinstance(col, bytes) else str(col) for col in tps_data_df.columns]

    # Set index to native datetime objects (NOT a DatetimeIndex)
    tps_data_df['BufTimes'] = buf_times
    tps_data_df.index = pd.Series(buf_times, dtype='object')  # This prevents automatic conversion to DatetimeIndex
    tps_data_df.drop(columns='BufTimes', inplace=True)

    return tps_data_df


#%% main functions

#%%load ltof data

def load_ltof_data(input_file_path_ltof, compounds=None, all_compounds=True):
    """
    Loads and concatenates data from multiple HDF5 files in the same folder.
    Ensures that all timestamps are preserved, filling missing values with NaN.
    Handles duplicate index values to avoid 'InvalidIndexError'.
    """
    # Create a list to store DataFrames
    dataframes = []

    # Iterate over all files in the LTOF data folder
    for file_name in os.listdir(input_file_path_ltof):
        # Construct the full path
        input_file_full_name = os.path.join(input_file_path_ltof, file_name)

        # Load the DataFrame from the HDF5 file
        peak_data_active_compounds_df = load_hdf_file(input_file_full_name, compounds, all_compounds)

        # Ensure index is in datetime format
        if not isinstance(peak_data_active_compounds_df.index, pd.DatetimeIndex):
            peak_data_active_compounds_df.index = pd.to_datetime(peak_data_active_compounds_df.index)

        # Remove any duplicate timestamps before resetting the index
        peak_data_active_compounds_df = peak_data_active_compounds_df[~peak_data_active_compounds_df.index.duplicated(keep="first")]

        # Reset index to allow safe concatenation
        peak_data_active_compounds_df = peak_data_active_compounds_df.reset_index()

        # Append to list
        dataframes.append(peak_data_active_compounds_df)

    # Concatenate all DataFrames, preserving all timestamps
    if dataframes:
        peak_data_df = pd.concat(dataframes, axis=0, join='outer', ignore_index=True)

        # Ensure index column exists
        if "index" in peak_data_df.columns:
            # Set back the DatetimeIndex after concatenation
            peak_data_df = peak_data_df.set_index("index")

        # Drop any duplicate timestamps again after merging
        peak_data_df = peak_data_df[~peak_data_df.index.duplicated(keep="first")]

        # Sort the index to maintain chronological order
        peak_data_df = peak_data_df.sort_index()
    else:
        peak_data_df = pd.DataFrame()  # Return an empty DataFrame if no data is loaded

    return peak_data_df




#%%load housekeeping data

#loads and concatenates data from several files in the same folder
def load_housekeeping_data(input_file_path_house_keeping):
    # Create an empty DataFrame to store the data     
    house_keeping_df = pd.DataFrame()
    #iterate over all files
    for file_name in os.listdir(input_file_path_house_keeping):
        # Load the text file into a DataFrame
        selected_data = pd.read_csv(os.path.join(input_file_path_house_keeping, file_name), sep=' ',index_col=0)
        house_keeping_df = pd.concat([house_keeping_df, selected_data])

    #Convert Mac HFS+ timestamp to datetime via unix conversion
    house_keeping_df.index = pd.to_datetime((house_keeping_df.index-2082844800),unit = 's')
    house_keeping_df['Housekeeping time'] = house_keeping_df.index
    
    return house_keeping_df

#%%Load MetNav data

def clean_file_remove_tabs(original_path, cleaned_path):
    with open(original_path, 'r', encoding='utf-8') as infile, open(cleaned_path, 'w', encoding='utf-8') as outfile:
        for line in infile:
            outfile.write(line.replace('\t', ''))

#highly custonized import due to individual data structure

def load_met_nav_data(input_file_path_met_nav, Metric=True):

    # Specified prefixes
    columns_to_keep = ['Time_Start', 'Latitude', 'Longitude', 'Pressure_Altitude', 'Static_Air_Temp', 'WindSpeed', 'WindDir', 'Wind_Speed', 'Wind_Dir']

    # Function to trim column names and rename specific columns
    def rename_column(col_name):
        # Rename "WindSpeed" and "WindDir" directly
        if col_name == "WindSpeed_m_s":
            return "Wind_Speed"
        elif col_name == "WindDir_deg":
            return "Wind_Dir"
        # For other columns, trim to the prefix
        for prefix in columns_to_keep:
            if col_name.startswith(prefix):
                return prefix  # Return the prefix as the new column name
        return col_name  # Return the original name if no prefix matches
    
    # Step 1: Clean file and remove all tabs
    cleaned_path = input_file_path_met_nav + ".cleaned"
    clean_file_remove_tabs(input_file_path_met_nav, cleaned_path)
    
    # Read the .ict file into a pandas DataFrame
    met_nav_df = pd.read_csv(cleaned_path, delimiter='\t')

    # Access the date and convert to datetime object
    entry_str = met_nav_df.iloc[5, 0]
    # Split the string based on commas
    parts = entry_str.split(',')
    # Extract year, month, and day
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    # Convert to datetime format
    date_time_obj = datetime(year, month, day)

    # Drop everything before the actuall data
    # Find indicies of all occurrences that start with 'Time_Start' 
    #(The second should be the line of the header no matter how much information was stored before)
    all_occurrences = met_nav_df[met_nav_df.iloc[:, 0].str.startswith('Time_Start')].index
    #the last occurance should eb the header of the table
    last_occurrence_index = all_occurrences[-1]
    #Drop all rows before the second occurrence
    met_nav_df = met_nav_df.loc[last_occurrence_index:]
    # Split into separate columns by seperating between commas
    met_nav_df = met_nav_df.iloc[:, 0].str.split(',', expand=True)
    # Move the first row values as column names
    met_nav_df.columns = met_nav_df.iloc[0]
    # Drop the first row
    met_nav_df.drop(met_nav_df.index[0], inplace=True)

    # Remove spaces from column labels
    met_nav_df.rename(columns=lambda x: x.replace(' ', ''), inplace=True)
    # Convert all entries to float format
    met_nav_df = met_nav_df.astype(float)
    # convert the 'Time_Start' columninto datetime fromat
    met_nav_df['Time_Start'] = date_time_obj + met_nav_df['Time_Start'].apply(lambda x: timedelta(seconds=x))

    #finter out only TS that are in all datasets
    filtered_columns = [col for col in met_nav_df.columns if any(col.startswith(prefix) for prefix in columns_to_keep)]
    # Creating a new DataFrame with only the filtered columns
    met_nav_df = met_nav_df[filtered_columns]
    #Apply renaming logic to the column names
    met_nav_df.columns = [rename_column(col) for col in met_nav_df.columns]

    #concert pressure altitude into meter if selected
    if Metric:
        met_nav_df['Pressure_Altitude'] = met_nav_df['Pressure_Altitude']*0.3048
        
    # Set the 'Time_Start' column as the index
    met_nav_df.set_index('Time_Start', inplace=True)
    #make a copy of the time for checking after merging
    met_nav_df['MetNav time'] = met_nav_df.index
    
    return met_nav_df


def load_env_flag_data(input_file_path_env_flag):

    # Specified prefixes
    columns_to_keep = ['Time_Start', 'leg_number', 'profile_number', 'BB_key', 'boundary_layer_key', 'urban_plume_key']

    # Read the .ict file into a pandas DataFrame
    env_flag_df = pd.read_csv(input_file_path_env_flag, delimiter='\t')

    # Access the date and convert to datetime object
    entry_str = env_flag_df.iloc[5, 0]
    # Split the string based on commas
    parts = entry_str.split(',')
    # Extract year, month, and day
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])
    # Convert to datetime format
    date_time_obj = datetime(year, month, day)
    
    #Drop all rows before the data starts
    env_flag_df = env_flag_df.loc[36:]
    # Split into separate columns by seperating between commas
    env_flag_df = env_flag_df.iloc[:, 0].str.split(',', expand=True)
    # Move the first row values as column names
    env_flag_df.columns = env_flag_df.iloc[0]
    # Drop the first row
    env_flag_df.drop(env_flag_df.index[0], inplace=True)

    # Remove spaces from column labels
    env_flag_df.rename(columns=lambda x: x.replace(' ', ''), inplace=True)
    # Convert all entries to float format
    env_flag_df = env_flag_df.astype(float)
    # convert the 'Time_Start' columninto datetime fromat
    env_flag_df['Time_Start'] = date_time_obj + env_flag_df['Time_Start'].apply(lambda x: timedelta(seconds=x))

    #finter out only TS that are in all datasets
    filtered_columns = [col for col in env_flag_df.columns if any(col.startswith(prefix) for prefix in columns_to_keep)]
    # Creating a new DataFrame with only the filtered columns
    env_flag_df = env_flag_df[filtered_columns]

    # Set the 'Time_Start' column as the index
    env_flag_df.set_index('Time_Start', inplace=True)
    #make a copy of the time for checking after merging
    env_flag_df['EnvFlag time'] = env_flag_df.index
    
    return env_flag_df
