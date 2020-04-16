import pymzml
import matplotlib.pyplot as plt
import numpy as np

import pyteomics
from pyteomics import mzml, auxiliary

import plotly.graph_objects as go

import pandas as pd

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

def motif_seek(ms2_scans, motifs, error = 0.002, noise_level = 10, precursor_base = 500, top_frags = 5, precursor_dist = 50, mzrange = [0, 500], rtrange = [0, 20]):
    motif_result = []
    motif_range = []
    for motif in motifs:
        motif_range.append([motif - error, motif + error])
    
    for scan in ms2_scans:
        precursor = scan.selected_precursors[0]['mz']
        drop_index = np.argwhere(scan.i <= noise_level)
        scan.i = np.delete(scan.i, drop_index)
        scan.mz = np.delete(scan.mz, drop_index)
        
        frag = scan.mz[scan.mz < precursor]
        frag_i = scan.i[: len(frag)] # In case need it
        base_index = np.argwhere(frag_i >= precursor_base)
        top_list = sorted(frag[base_index])[-top_frags : ]
        top_range = precursor - precursor_dist
        top_list = [top for top in top_list if top >= top_range]
        
        
        neutral_loss = precursor - frag
        for top_frag in top_list:
            neutral_loss_top = top_frag - frag
            neutral_loss_top = neutral_loss_top[neutral_loss_top > 0]
            neutral_loss = np.append(neutral_loss, neutral_loss_top)
        
        mtf_count = 0
        for mtf in motif_range:
            mtf_hit = neutral_loss[(mtf[0] < neutral_loss) & (neutral_loss < mtf[1])]
            if mtf[0] == motifs[0] - error: #Only show the mtf_error for first motif
                mtf_error = mtf_hit - motifs[0]
            if len(mtf_hit) > 0:
                mtf_count += 1
        
        if mtf_count == len(motifs):
            motif_result.append([scan.selected_precursors[0]['mz'], round(scan.scan_time[0],2), scan.ID, np.round(mtf_error, 4)])
    
    result = sorted(motif_result)
    d_motif = pd.DataFrame(result, columns = ['mz', 'rt', 'scanID','motif_error']) 
    d_motif = d_motif[(d_motif['mz'] < mzrange[1]) 
                    & (d_motif['mz'] > mzrange[0])
                    & (d_motif['rt'] > rtrange[0])
                    & (d_motif['rt'] < rtrange[1])]
    d_motif.sort_values(by=['mz', 'rt'], inplace=True)
    d_motif.reset_index(inplace=True)
    d_motif.drop(columns = ['index'], inplace=True)

    return d_motif

def motif_export(motif_result, export_name='Motif.xlsx'):
    writer = pd.ExcelWriter(export_name, engine='xlsxwriter')
    for i,j in motif_result.groupby('mz'):
        d_iter = pd.DataFrame(j)
        d_iter.to_excel(writer, sheet_name=str(round(i,2)))
    writer.save()
    
    return

def find_scan(ms2s, scanid, interactive = True):
    
    for scan in ms2s:
        if scan.ID == scanid:
            break
    
    print('Precursor m/z: {:0.2f}, Scan time: {:0.1f} minute'.format(scan.selected_precursors[0]['mz'], scan.scan_time[0]))
    
    mz = scan.mz
    ints = scan.i
    
    if interactive == True:
        plt.clf()
        fig = go.Figure([go.Bar(x=mz, y=ints, marker_color = 'red', width = 0.5,
                        hovertemplate =
                        'Int: %{y}'+
                        '<br>m/z: %{x}<br>')])
        fig.update_traces(marker_color='rgb(158,202,225)', marker_line_color='rgb(0,0,0)',
                  marker_line_width=0.5, opacity=1)
        fig.update_layout(
                template = 'simple_white',
                width = 1000,
                height = 600,
                xaxis = dict(title = 'm/z ratio',
                        rangeslider=dict(
            visible = True
        )),
                yaxis = dict(
                    title = 'Intensity'))
        fig.show()
    
    elif interactive == False:
        plt.figure(figsize=(10,5))
        plt.bar(mz, ints, width = 1.0)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        plt.title('MS1 spectrum')
        plt.xlim(0,340)

    return 

def batch_scans(path, remove_noise = True, thres_noise = 1000):
    all_files = glob.glob(path + "/*.mzML")
    scans = []
    file_list = []
    for file in tqdm(all_files):
        scan = get_scans(file)
        if remove_noise == True:
            noise_removal(scan, thres_noise)
        scans.append(scan)
        file_list.append(Path(file).name)
    print(file_list)
    print('Batch read finished!')
    
    return scans, file_list