import pymzml
import numpy as np

import pandas as pd
import glob
from tqdm import tqdm
from pathlib import Path 

from datetime import datetime
now = datetime.now()
dt_string = now.strftime("%d%m%Y %H-%M-%S")

def get_scans(path, ms_all = False, ms_lv = 1):
    run = pymzml.run.Reader(path)
    scans = []
    if ms_all == False:
        for scan in run:
            if scan.ms_level == ms_lv:
                scans.append(scan)
    elif ms_all == True:
        for scan in run:
            scans.append(scan)
            
    return scans

def ms_qc(mzml_scans):
    # seperate ms1 and ms2
    ms1 = []
    ms2 = []
    for scan in mzml_scans:
        if scan.ms_level == 1:
            ms1.append(scan)
        elif scan.ms_level == 2:
            ms2.append(scan)
    
    # ms1 QC
    tot_ms1 = len(ms1)
    TIC = []
    for scan in ms1:
        TIC.append(scan.TIC)
        
    tic_h = max(TIC)
    tic_avg = sum(TIC) / len(TIC)
    
    #ms2 QC
    if len(ms2) != 0:
        tot_ms2 = len(ms2)
        h_response = []
        precursor_mz = []
        precursor_i = []
        for scan in ms2:
            h_response.append(scan.i.max())
            precursor_mz.append(scan.selected_precursors[0]['mz'])
            precursor_i.append(scan.selected_precursors[0]['i'])
        avg_response_h = sum(h_response) / len(h_response)
        avg_precursor_i = sum(precursor_i) / len(precursor_i)
        l_precursor_i = min(precursor_i)
        precursor_range = str(round(min(precursor_mz), 2)) + '~' + str(round(max(precursor_mz),2))
    
    result = [tot_ms1, "{:.2e}".format(tic_h), "{:.2e}".format(tic_avg), tot_ms2, round(avg_response_h, 2), "{:.2e}".format(avg_precursor_i), round(l_precursor_i, 2), precursor_range]

    return result

def batch_scans(path):
    all_files = glob.glob(path + "/*.mzML")
    scans = []
    file_list = []
    for file in tqdm(all_files):
        scan = get_scans(file, True)
        scans.append(scan)
        file_list.append(Path(file).name)
    print(file_list)
    print('Batch read finished!')
    
    return scans, file_list

def qc_gen(path):
    batch_scan, file_list = batch_scans(path)
    print('All files read in!')
    
    qc_result = []
    ran_scan = len(batch_scan)
    print('Generating QC report...')
    for index in tqdm(np.arange(ran_scan)):
        result = ms_qc(batch_scan[index])
        
        result = [file_list[index]] + result
        qc_result.append(result)
    print('Generating dataframe...')
    col = ['file_name', 'total_ms1_scan', 'max_tic', 'avg_tic','tot_ms2', 'avg_ms2_max', 'avg_ms2precursor_i', 'min_ms2precursor_i', 'ms2precursor_range']
    d_result = pd.DataFrame(qc_result, columns = col)
    print('Generating report to csv...')
    
    return d_result

mzml_path = input('Please input mzml path:')
d_result = qc_gen(mzml_path)
name = dt_string + 'QC report.csv'
d_result.to_csv(name)
print('Finished!')