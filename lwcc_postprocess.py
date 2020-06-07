#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculation and salinity correction of absorption spectra by CDOM from 
LWCC .raw files.

                                                                     
About LWCC sample labeling: 
Sample labels before the separator (e.g. '#')) need to contain certain keywords:
underway sample - keyword "UW" (case insensitive),
NaCl solution - keyword "NaCl", "PSU" or "Sal" (case insensitive),
and CTD sample labels do not contain any keywords mentioned in the other two catalogs.
                                                                     

Before using this script, please 
1) install:
- Python 3.7x or higher version    
- Essential Python packages: numpy, scipy, pandas, matplotlib.
2) prepare the following files (see examples as well as the main program) and 
put them in the same directory as the one containing .raw files:
- config_PostProc.txt
- rmSpectraIndex.txt
- label_underway.txt and/or label_station.txt
- lwcc_labels_modify.txt (no need if labelling was done as mentioned above.)
- labels_latlon_datetime_lwcc.txt (no need if data are not meant for uploading to Pangaea.)

Analysis method detailed in:
1) Lefering et al. (2017). Uncertainty budgets for liquid waveguide CDOM 
absorption measurements. Applied Optics, 56(22), 6357-6366.
2) Röttgers (2020). In prep.
3) Dall’Olmo et al. (2017). Determination of the absorption coefficient of 
chromophoric dissolved organic matter from underway spectrophotometry. Optics 
Express, 25(24), A1079-A1095.

@author: 
    Rüdiger Röttgers (ruediger.roettgers@hzg.de) (raw data process);
    Yangyang Liu (yangyang.liu@awi.de) (postprocess), February 2020.
"""

import glob, os, shutil, math, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as inp
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def readConfig(config):
    paramDict = {}
    try:
        with open(config) as f:
            for line in f:
                if line.startswith('#')==False and len(line)>1:
                    key = line.strip().split('=')[0]
                    try:
                        value = line.strip().split('=')[1].strip()
                    except:
                        value = None
                    
                    paramDict[key] = value
    except Exception as e:
        print(e)
    finally:
        return paramDict

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------   
def extract_pangaea_tempsal(tsg_pangaea, station=True):
    '''
    This function extracts temperature and salinity data downloaded from Pangae. 
    
    Input:
    1) tsg_pangaea - str. Path of the TS data file "*oce.tab".
    2) station - bool, default True. If True, profile TS data; if False, 
    underway TS data.
    
    Output:
    Extracted TS data saved as "*oce_extracted.txt" (tab delimited), 
    in which the first row is header. For profile data, the columns are label, 
    datetime, depth, temperature and salinity, respectively. For underway data, 
    the columns are datetime, temperature and salinity, respectively.
    '''
    foo = tsg_pangaea.split('/')
    print(f'Importing {foo[-1]}')

    with open(tsg_pangaea, 'r') as f:
        for num, line in enumerate(f,1):
            if '*/' in line:
                break  
        columns = f.readlines()[0]
        
    columns = columns.split('\t')
    data = pd.read_csv(tsg_pangaea, names=columns, skiprows=num+1, sep='\t') 

    
    if station: #CTD station
        data_new = pd.DataFrame()
        data_new['label'] = data[[s for s in columns if "Event" in s][0]]
        data_new['datetime'] = data[[s for s in columns if "Date" in s][0]]
        data_new['depth'] = data[[s for s in columns if "Depth" in s][0]]
        data_new['temperature'] = data[[s for s in columns if "Temp" in s][0]]
        data_new['salinity'] = data[[s for s in columns if "Sal" in s][0]]
        data_new['depth'] = round(data_new['depth'])
    else: #underway
        data_new = data.drop(['Latitude', 'Longitude', 'Depth water [m]'], axis=1)
        data_new.columns=['datetime', 'temperature', 'salinity']
    
    data_new = data_new.interpolate(method='nearest')       
    data_new.to_csv(tsg_pangaea.replace('.tab','_extracted.txt'), 
                        index=False, header=True, encoding='utf-8', sep='\t')
    print(f'Temperature and salinity data from {tsg_pangaea} Extracted!') 
    
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------- 
class lwccPostProc():
    
    def __init__(self, config='config_PostProc.txt'):
        paramDict = readConfig(config)
        xi = paramDict['xi']
        wlstart, increment, wlend = map(float, xi.split(':'))
        self.xi = np.arange(wlstart, wlend+increment, increment)
        self.path_length=float(paramDict['path_length'])
        self.separator = paramDict['separator']
        self.tsgST_pangaea = paramDict['tsgST_pangaea']
        self.tsgUW_pangaea = paramDict['tsgUW_pangaea']
        self.labelST = paramDict['labelST']
        self.labelUW = paramDict['labelUW']
        self.TFpangaea = eval(paramDict['TFpangaea'])
        self.concurrent_info = paramDict['concurrent_info']
        self.modified_labels = paramDict['modified_labels']
    
    #--------------------------------------------------------------------------
    def mk_spec(self, line, line_dark):
        """ 
        This function subtracts a dark measurement from either a reference or sample 
        measurement, and interpolate the data in the range of "xi".
        """
    
        data = line.split('[ms]=')[1].split()
    #    inttime = float(data[0])
        spec = np.array(data[1:], dtype='float')
        
        data_dark = line_dark.split('[ms]=')[1].split()
    #    inttime = float(data_dark[0])    
        spec_dark = np.array(data_dark[1:], dtype='float')
        
        spec = spec - spec_dark
        spec = spec - np.mean(spec[0:5]) #simple straylight correction, subtraction of constant counts
      
        if len(self.xi) == len(spec):
            return spec
        elif  len(self.wl_raw) == len(spec):
            #print ('#interpolate')
            f = inp.interp1d(self.wl_raw, spec)
            return f(self.xi)
        else:
            print ('>>mk_spec ERROR! Lengths of measured spectrum and xi differ!')
        
    #--------------------------------------------------------------------------
    def output(self, fw, line_spec, spectra):
        '''
        This function writes spectral data in a file fw.
        '''
        
        labels = line_spec.split()
        fw.write('%s %s %s'%tuple(labels[:3]))
        for spec in spectra:
            fw.write(' %.6f'%spec)
        fw.write('\n')
    
    #--------------------------------------------------------------------------
    def processRaw(self, filename):
        base = os.path.basename(filename)  
        print('Processing ', base)
        filename_l0 = base.replace('.raw','.l0a')   
        with open(filename, 'r') as f:
            lines = f.readlines()
        with open(filename_l0, 'w') as fw: 
            for i in range(len(lines)):
                line = lines[i]
                if 'Wavelen'in line:
                    self.wl_raw = np.array(line.split(':')[1].split(), dtype='float')
                    fw.write('%'+'Pathlength used for calculations: %2.4f m'%self.path_length)
                    fw.write('\n')
                    fw.write('%wavelength:')
                    for wl in self.xi: 
                        fw.write(' %i'%wl)
                    fw.write('\n')    
                elif line.startswith('%'):
                    fw.write(line)
                elif 'dark' in line:   
                    idark = i
                    if i >= len(lines): break
                elif 'ref' in line: 
                    iref = i
                    try:
                        ref_spec = self.mk_spec(lines[iref], lines[idark])
                    except:
                        print('Warning: data should start with dark measurements! \n')   
                    if i >= len(lines): break
                elif 'spec' in line:
                    isample = i
                    try:
                        sam_spec = self.mk_spec(lines[isample], lines[idark])
                        #calcualte absorption from light intensity
                        #np.log - natural logarithm
                        absor = -np.log(sam_spec/ref_spec)/self.path_length
                        self.output(fw, line, absor)
                    except:
                        pass
                        print('Warning: a sample spectrum may have WRONG \
                              length and was removed! Check if there is any \
                              data loss in .l0a file!\n')
                    if i >= len(lines): break
                else:
                    print('#no match', i+1, line[:80])
        
    #--------------------------------------------------------------------------
    def getWL(self, filename):
        with open(filename, 'r') as f:
            for line in f:
                if 'wavelength' in line:
                   break        
        wl = pd.Series(line.split(' ')[1:]).astype('float')
        
        wl_selected = [300, 350, 400, 430, 440, 490, 650, 690, 700, 710, 720, 
                       730, 750]
        wlpos = {}
        for i, wvl in enumerate(wl_selected):
            key = 'pos'+str(wvl)
            value = np.where(wl>=wvl)[0][0]
            wlpos[key] = value
            
        return wl, line, wlpos
    
    #--------------------------------------------------------------------------
    def extract_abs(self, filename):
        
        data = pd.read_csv(filename, header=None, comment='%', sep=' ') 
        label = data.iloc[:,2]
        label2 = label.copy()
        for i,lb in enumerate(label):
            if (self.separator !='') & (self.separator in lb):
                label2.iloc[i] = lb[:lb.find(self.separator)]
            else:
                label2.iloc[i] = lb
        label = label2 
        
        return data, label
    
    #--------------------------------------------------------------------------
    def plotSpectrum(self, dfx, dfy, figname, ylabel, xlabel='Wavelength [$nm$]', 
                     legend=None):
        fig = plt.figure(figsize=(8.5, 6.5))
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.tick_params(axis="both", labelsize=10)
        ax.plot(dfx, dfy, linewidth=1)
        if legend is not None:
            ax.legend(legend, loc='upper left', bbox_to_anchor=(1, 1),
                      fontsize='small')
        ax.set_xlabel(xlabel,fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)
        fig.savefig(figname, dpi=200)
        plt.close(fig)
        
    #--------------------------------------------------------------------------
    def saveData(self, filename, dfvarname, headline=None, index=True, 
                 header=False, sep='\t'):
        
        with open(filename, 'w') as f:
            if headline is not None:
                f.write(headline)
            dfvarname.to_csv(f, index=index, header=header, sep=sep, 
                             encoding='utf-8') 
    
    #--------------------------------------------------------------------------
    def createDir(self, dirname, overwrite=False):
        
        if overwrite:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            else:
                shutil.rmtree(dirname) #removes all the existing directories!
                os.makedirs(dirname)
        else:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
    
    #--------------------------------------------------------------------------
    def plot_l0a(self, filename):
        '''
        This function plots LWCC absorption data of each sample in a .l0a file.
        
        Input:
        filename - str. Path of a LWCC .l0a file.
        Output:
        A folder named as the base name of the .l0a file in the working direcoty 
        containing plots from individual sample measurements. The legend of each 
        absorption spectrum in every plot gives the row index of each data point 
        when this .l0a data file is imported into a pandas dataframe via 
        "pd.read_csv". These row indices can be used to prepare the file 
        "rmSpectraIndex.txt". Note: if the number of samples is greater than 4, 
        check the labelling in .raw data!
        '''

        foo = filename.split('/')
        print(f'Importing {foo[-1]}...')
        
        #get wavelengths
        wl,line, wlpos = self.getWL(filename)       
        try:
            #Extract absorption spectra
            data, label = self.extract_abs(filename)
            label_uni = np.unique(label,return_index=True, return_inverse=True, 
                               return_counts=True)
            #Create target Directory 
            dir_fig = foo[-1][:-4]
            self.createDir(dir_fig, overwrite=True)

            #plot absorption spectra
            for i,lb in enumerate(label_uni[0]):
                pos = np.where(label_uni[2] == i)[0]            
                data_for_plot = data.iloc[pos, 3:].transpose()
                lgd = [str(idx) for idx in pos]
                figname = os.path.join(dir_fig, 'check_'+ lb + '.png')
                self.plotSpectrum(wl, data_for_plot, figname, legend=lgd,
                                  ylabel='$Level\ 0\ uncorrected\ a_{CDOM}$') 
            print('Plotting Done!') 
        except:
            print(f'Warning: No data found in: {foo[-1]}!!! \n')

    #--------------------------------------------------------------------------    
    def rmabs_median(self, filename, row_index, plot=True):  
        '''
        This function removes suspicious measurements in LWCC .l0a data indexed 
        by "row_index" and take the median values and standard deviation of 
        repeatedly measured data.
        
        Input:
        1) filename - str. Path of a LWCC .l0a file.
        2) row_index - list, default None. List of numbers as row indices of 
        suspicious data points when this .l0a data file is imported into a pandas 
        dataframe via "pd.read_csv". The indices of each data point are annotated 
        in the legend of the plots produced by the function "plot_l0a".
        3) plot - bool, default True. Plot the output spectra.
        Output:
        1) Tab delimited files in the direcory "Median_agfiles" named as "*_ag_median.txt" 
        and "*_ag_sd.txt", respectively, with "*" standing for the base name of .l0a files.
        2) Plots of the spectra in the "*_ag_median.txt" and "*_ag_sd.txt" named as 
        "*_ag_median.png" and "*_ag_sd.png", respectively, in the directory produced 
        by the function "plot_l0a". 
        '''
        
        filename = filename.strip()
        foo = filename.split('/')
        dir_fig = foo[-1][:-4]
        print(f'Importing {foo[-1]}...')
        
        wl, line, wlpos = self.getWL(filename)
        #Extract absorption spectra data
        data, label = self.extract_abs(filename)
        data.iloc[:,2] = np.array(label)
        
        if not math.isnan(row_index[0]):
            data.loc[row_index, :] = np.nan
        
        #take the median of repeated measurements in each file
        data_median = data.groupby(pd.Grouper(key=2)).median() 
        data_median.iloc[:,0] = data_median.index
                
        #Create target Directory 
        self.createDir('Median_agfiles')
        
        outname_median = filename.replace('.l0a','_ag_median.txt') 
        self.saveData(os.path.join('Median_agfiles',outname_median), data_median, 
                      headline=line, index=False) 
       
        #take the standard deviation of repeated measurements in each file
        data_sd = data.groupby(pd.Grouper(key=2)).agg(np.std, ddof=1)
        data_sd.iloc[:,0] = data_sd.index
        data_sd = data_sd.drop(columns=[1])

        outname_sd = filename.replace('.l0a','_ag_sd.txt')  
        self.saveData(os.path.join('Median_agfiles',outname_sd), data_sd, 
                      headline=line, index=False) 
        
        if plot:
            figname = os.path.join(dir_fig, dir_fig + '_ag_median.png')
            #plot acdom spectra
            acdom = data_median.iloc[:,1:]
            wl_for_plot = wl.iloc[wlpos['pos300']:wlpos['pos720']+1]
            acdom_for_plot = acdom.iloc[:,wlpos['pos300']:wlpos['pos720']+1].transpose()
            lgd = data_median.iloc[:,0].to_list()
            if len(acdom) != 0:
                self.plotSpectrum(wl_for_plot, acdom_for_plot, figname, legend=lgd,
                                  ylabel='$Uncorrected\ a_{CDOM}$ [$m^{-1}$]')
            else:
                print(f'No plot: after removing suspicious data, no data were left \
                      in {foo}!')
        
        print(f'Removal of suspicious CDOM absorption coefficients and averaging for {foo[-1]} Finished! \n') 
    
    #--------------------------------------------------------------------------
    def get_labels(self):
        '''
        This function extracts sample labels from all '*_ag_median.txt' files.
        Its output can be used to construct the file for the parameter 
        "modified_labels" in the configuration file "config_PostProc.txt".
        Output:
        A file named as 'lwcc_merged_labels.txt' containing all sample labels.
        '''
        print ('Getting the list of LWCC labels...')
        fnames=sorted(glob.glob(os.path.join('Median_agfiles',
                                             '*ag_median.txt'))) 
        
        absorption = pd.DataFrame()
        for f in fnames:
            #Read data file
            print(f'Reading {f} \n')     
            try:
                data = pd.read_csv(f, header=None, comment='%', sep='\t')
                absorption = absorption.append(data)
                del data
            except:
                pass
                print(f'Warning: No data found in: {f}!!!')      
        
        #take the median of repeated measurements (if any)
        absorption = absorption.groupby(pd.Grouper(key=0)).median()   
        absorption.insert(0,"label",absorption.index, True) 
        absorption.sort_index(inplace=True)
        label = absorption.iloc[:,0]
        label.to_csv('lwcc_merged_labels.txt', index=False, header=False, 
                     encoding='utf-8')

    #--------------------------------------------------------------------------    
    def merge_abs(self, keyword='median', plot=True):   
        '''
        This function merges all '*_ag_median.txt' files, seperates data from CTD 
        sations, underway and NaCl solution, and merges data from each catalog, 
        respectively.
        
        The prerequisit for data seperation is that underway sample labels contain 
        the keyword "UW" (case insensitive), NaCl solution labels contain 
        keywords "NaCl", "PSU" or "Sal" (case insensitive), and CTD sample labels
        do not contain any keywords mentioned in the other two catalogs. If the 
        sample naming does not follow this, modify the labels by use of parameter 
        "modified_labels" in the configuration file "config_PostProc.txt".
        
        Input:
        1) keyword - str, default 'median'. 'meidan' and 'sd' stand for the files to be 
        merged are median or standard deviation of repeterd absorption spectra 
        measurements, respectively.
        2) plot - bool, default True. Plot all sample spectra.     
        Output:
        1) Eight tab delimited files saved in the working directory as 
        "lwcc_merged_median_ag_salt.txt" (absorption spectra by NaCl solution), 
        "lwcc_merged_sd_ag_salt.txt" (standard deviation of absorption spectra by NaCl solution), 
        "lwcc_merged_median_ag_underway.txt" (absorption spectra by underway samples), 
        "lwcc_merged_sd_ag_underway.txt" (standard deviation of absorption spectra by underway samples), 
        "lwcc_merged_median_ag_station.txt" (absorption spectra by CTD samples),
        "lwcc_merged_sd_ag_station.txt" (standard deviation of absorption spectra by CTD samples), and
        "lwcc_merged_median_ag_all.txt" (absorption spectra by all samples), 
        "lwcc_merged_sd_ag_all.txt" (standard deviation of absorption spectra by all samples), 
        respectively.
        2) Two plots saved in the working directory as "lwcc_merged_median_ag_all.png" 
        (all absorption spectra) and "lwcc_merged_sd_ag_all.png" (standard deviation 
        of all absorption spectra), respectively.
        '''
        print ('Merging CDOM absorption data files...')
        fnames = sorted(glob.glob(os.path.join('Median_agfiles',
                                             '*ag_'+keyword+'.txt'))) 
        wl, line, wlpos = self.getWL(fnames[0])    
        absorption = pd.DataFrame()
        for file in fnames:
            #Read data file
            print(f'Reading {file} \n')     
            try:
                data = pd.read_csv(file, header=None, comment='%', sep='\t')
                absorption = absorption.append(data)
                del data
            except:
                pass
                print(f'Warning: No data found in: {file}!!!')   
        
        if self.modified_labels != '':
            new_labels = pd.read_csv(self.modified_labels, sep='\t')
            for i, lb_ori in enumerate(absorption.iloc[:,0]):
                for j, lb_new in enumerate(new_labels.iloc[:,0]):
                    if lb_new == lb_ori:
                       absorption.iloc[i,0] = new_labels.iloc[j,1]     
        
        #uniformly relabel NaCl as 'NaCl_91psu'
        absorption.index = range(len(absorption))
        label = absorption.iloc[:,0]
        idx_salt = label.str.contains('PSU|Sal|NACL|NACI', case=False)
        idx_salt = list(idx_salt[idx_salt==True].index)
        absorption.iloc[idx_salt,0] = 'NaCl_91psu'
        
        #take the median of repeated measurements (if any)
        absorption = absorption.groupby(pd.Grouper(key=0)).median() 
        absorption.sort_index(inplace=True)
        absorption.insert(0,"label",absorption.index, True) 
        absorption.index = range(len(absorption))
        label = absorption.iloc[:,0]
        
        #absorption of NACL solution 
        absorption_salt = absorption[label.str.contains('PSU|Sal|NACL|NACI', 
                                                        case=False)]
        #absorption of all samples
        absorption = absorption[~absorption['label'].str.contains(
                'PSU|Sal|NACL|NACI|MQ|SW', case=False)]
            
        #absorption of underway samples
        absorption_uw = absorption[label.str.contains('UW', case=False)]

        #absorption of CTD station samples
        absorption_station = absorption[~absorption['label'].str.contains('PSU|Sal|NACL|NACI|UW|MQ|SW',
                                    case=False)]
    
        outnames = ['lwcc_merged_'+keyword+'_ag_salt_uncorrected.txt', 
                'lwcc_merged_'+keyword+'_ag_underway_uncorrected.txt',
               'lwcc_merged_'+keyword+'_ag_station_uncorrected.txt', 
               'lwcc_merged_'+keyword+'_ag_all_uncorrected.txt']
        varnames = ['absorption_salt', 'absorption_uw', 'absorption_station', 
                'absorption']
        
        for i in range(len(outnames)):  
            varname = eval(varnames[i])
            self.saveData(outnames[i], varname, headline=line, index=False)
        if plot:
            acdom = absorption.iloc[:,1:]
            wl_for_plot = wl.iloc[wlpos['pos300']:wlpos['pos750']+1]
            acdom_for_plot = acdom.iloc[:,wlpos['pos300']:wlpos['pos750']+1].transpose()
            figname = 'lwcc_merged_'+keyword+'_ag_all.png'
            self.plotSpectrum(wl_for_plot, acdom_for_plot, figname, 
                              ylabel='$Uncorrected\ a_{CDOM}$ [$m^{-1}$]')
            
        print ('Merging CDOM absorption data files finished!')

    #--------------------------------------------------------------------------        
    def abs_salcorr0(self, filename, filename_tsg, filename_label, 
                     smooth=True, station=True, plot=True):  
        '''
        This function corrects salinity effects and interference pattern by smoothing
        and offset in the NIR for LWCC CDOM absorption data (Lefering et al., 2017; 
        Röttgers (2019), in prep). When a spectrum is significantly positive in 
        the NIR (e.g. > 0.05 m^-1 in this case), it is not considered as a valid 
        measurement and removed.
        
        Input:
        1) filename - str. Path of CDOM absorption data file (tab delimited).
        2) filename_tsg - str. Path of Temperature and salinity data file (tab 
        delimited). Inside, the first row is header; the columns for profile data 
        are label, datetime, depth, temperature and salinity, respectively, and 
        for underway data are datetime, temperature and salinity, respectively.
        3) filename_label - str. Path of a tab delimited .txt file in which the 
        columns for profile data are LWCC labels, standard station labels in TSG 
        file and depth in meter, respectively, and for underway data are LWCC 
        labels, date (yyyy-mm-dd) and time (HH:MM), respectively.
        4) smooth - bool, default True. Method: Savitzky Golay first order 
        polynomial filter (window size, 5).
        5) station - bool, default True. If True, profile data; if False, underway 
        data.
        6) plot - bool, default True. Plot salinity corrected absorption spectra.
        Output:
        1) Salinity effect correcteted LWCC CDOM absorption data saved as 
        "*salcorr_smooth_null_final.txt". 
        2) Salinity effect correcteted LWCC CDOM absorption data plotted as
        "*salcorr_smooth_null_final.png".
        3) Individual salinity effect correcteted spectra plotted and saved in
        the directory "*salcorr_smooth_null".
        '''
        
        #Salinity correction coefficient
        filename_nacl = 'lwcc_merged_median_ag_salt_uncorrected.txt'
        data_nacl = pd.read_csv(filename_nacl, header=None, comment='%', sep='\t')
        #plt.plot(wl, data_nacl.squeeze().iloc[1:])
        
        #NaCl solution was prepared by adding 100 mL purified water to 10g NaCl,
        #and the salinity of the NaCl solution is about 91 PSU (g∕kg). 
        psi_sal =  data_nacl.iloc[:,1:] / 91 
        
        wl, line, wlpos = self.getWL(filename)
        #Salinity correction coefficient
        data_nacl = pd.read_csv(filename_nacl, header=None, comment='%', sep='\t')
        #plt.plot(wl, data_nacl.squeeze().iloc[1:])
        
        #NaCl solution was prepared by adding 100 mL purified water to 10g NaCl,
        #and the salinity of the NaCl solution is about 91 PSU (g∕kg). 
        psi_sal =  data_nacl.iloc[:,1:] / 91 
        
        data = pd.read_csv(filename, header=None, comment='%', sep='\t')
    #    label = data.iloc[:,0]
    #    depth = [label[i].split('#')[-1][:-1] for i in range(len(label))]


        #match TSG salinity data with LWCC data.
        data_tsg = pd.read_csv(filename_tsg, header=None, sep='\t', skiprows=1)

        if station: #CTD station
            
            figname = 'lwcc_merged_median_ag_station_salcorr.png'
        
            label_tsg = data_tsg[0].to_list() 
            label_tsg_uni = np.unique(label_tsg) 
            #cut the end of TSG labels, e.g. "-1", "-3", "-8"       
            label_tsg_uni_cut = [label_tsg_uni[i][:label_tsg_uni[i].find('-')] for i 
                                               in range(len(label_tsg_uni))] 
            label_tsg_uni_cut_uni = np.unique(label_tsg_uni_cut) 
            
         
            pos_same = []
            pos = []
            for lb in label_tsg_uni_cut_uni:
                for i in range(len(label_tsg_uni_cut)):
                    if label_tsg_uni_cut[i] == lb:
                        pos_same.append(i)
                #if there are same cutted unique labels, e.g. if there are two
                #"PS107_12" labels, take the last one.
                pos.append(pos_same[-1])
            label_tsg_uni_final = label_tsg_uni[pos]
            
            label_tsg_uni = np.vstack([label_tsg_uni, label_tsg_uni_cut]).transpose()
            label_tsg_uni_final = np.vstack([label_tsg_uni_final, label_tsg_uni_cut_uni]).transpose()
            del label_tsg_uni_cut, label_tsg_uni_cut_uni, pos_same, pos, i, lb
    
            stationlabel = pd.read_csv(filename_label, sep='\t') 
            
            label_station = stationlabel.iloc[:,1].to_list()
            label_station_cut = [label_station[i][:label_station[i].find(
                    '-')] for i in range(len(label_station))]
        
            label_station_uni, indices, counts = np.unique(label_station, return_index=True, 
                                                   return_counts=True)    
            label_station_uni_cut = [label_station_uni[i][:label_station_uni[i].find(
                    '-')] for i in range(len(label_station_uni))]
            label_station_uni = np.vstack([label_station_uni, label_station_uni_cut]).transpose()
            label_station = np.vstack([label_station, label_station_cut]).transpose()
            del label_station_uni_cut, label_station_cut
            
            ts_match = pd.DataFrame()
            for i in range(len(label_station_uni)):
                lb = label_station_uni[i,1]
                for j in range(len(label_tsg_uni_final)):
                    if label_tsg_uni_final[j,1]==lb:
                        label_ori = str(label_tsg_uni_final[j,0])
                        pos = [j for j,lbts in enumerate(label_tsg) if lbts==label_ori]
                        data_tsg_tmp = data_tsg.iloc[pos,:]
                        depth_station = stationlabel.iloc[indices[i]:indices[i]+
                                                          counts[i],:].iloc[:,2]
                        depth_station = round(depth_station)
                        #interpolate TSG depth into 1-meter interval
                        depth_tsg_interp = np.arange(min(min(data_tsg_tmp.iloc[:,2]),min(depth_station)), 
                                                     max(max(depth_station),max(data_tsg_tmp.iloc[:,2]))+1, 1)
                        
                        func_temp = inp.interp1d(data_tsg_tmp.iloc[:,2], data_tsg_tmp.iloc[:,3], 
                                            'nearest', fill_value='extrapolate')
                        func_sal = inp.interp1d(data_tsg_tmp.iloc[:,2], data_tsg_tmp.iloc[:,4], 
                                            'nearest', fill_value='extrapolate')
                        ts_interp = np.vstack([depth_tsg_interp, func_temp(depth_tsg_interp),
                                               func_sal(depth_tsg_interp)]).transpose()
                        
                        ts_index = [np.where(depth_tsg_interp==dep)[0][0] for dep in depth_station]
                        ts_match0 = pd.DataFrame(ts_interp[ts_index], index=depth_station.index)
                        ts_match = ts_match.append(ts_match0)
                        del depth_tsg_interp, func_temp, func_sal, ts_interp, label_ori
                        del depth_station, ts_index, ts_match0
                del lb
        
            stationlabel = pd.concat([stationlabel, ts_match], axis=1)
            stationlabel.index = stationlabel.iloc[:,0]
            stationlabel.iloc[:,1] = stationlabel.iloc[:,2]
            label_station = list(stationlabel.iloc[:,0])
    
        
        else: #underway
            
             figname = 'lwcc_merged_median_ag_underway_salcorr.png'
             stationlabel = pd.read_csv(filename_label, sep='\t', 
                                        parse_dates=[[1, 2]]) 
             
             #interpolate TSG time into 1-minute interval
             data_tsg.iloc[:,0] = pd.to_datetime(data_tsg.iloc[:,0])
             
             if min(stationlabel.iloc[:,0]) < min(data_tsg.iloc[:,0]):
                 pos_min = np.where(data_tsg.iloc[:,0] == min(data_tsg.iloc[:,0]))[0]
                 data_tsg = data_tsg.append({0:min(stationlabel.iloc[:,0]), 
                                  1:data_tsg.values[pos_min,1][0], 
                                  2:data_tsg.values[pos_min,2][0]}, ignore_index=True)
             
             if max(stationlabel.iloc[:,0]) > max(data_tsg.iloc[:,0]):
                 pos_max = np.where(data_tsg.iloc[:,0] == max(data_tsg.iloc[:,0]))[0]
                 data_tsg = data_tsg.append({0:max(stationlabel.iloc[:,0]), 
                              1:data_tsg.values[pos_max,1][0], 
                              2:data_tsg.values[pos_max,2][0]}, ignore_index=True)
                 
             data_tsg.index = data_tsg.iloc[:,0]    
             data_tsg = data_tsg.sort_index()
                 
             data_tsg_upsampled = data_tsg.resample('1min')
             data_tsg_interp = data_tsg_upsampled.interpolate(method='nearest') 
             data_tsg_interp.iloc[:,0] = data_tsg_interp.index
    
             #match station time with TSG time
             time_station = list(stationlabel.iloc[:,0])
             time_tsg = list(pd.to_datetime(data_tsg_interp.iloc[:,0]))
             time_overlap = list(set(time_station).intersection(time_tsg))
             idx_time_station = [time_station.index(dt) for dt in time_overlap]
             idx_time_tsg = [time_tsg.index(dt) for dt in time_overlap]
             stationlabel_dtoverlap = stationlabel.iloc[idx_time_station,:]
             data_tsg_dtoverlap = data_tsg_interp.iloc[idx_time_tsg,:]
             stationlabel_dtoverlap.index = data_tsg_dtoverlap.index
             stationlabel = pd.concat([stationlabel_dtoverlap, 
                                       data_tsg_dtoverlap], axis=1, sort=True)
             
             label_station = list(stationlabel.iloc[:,1])
             stationlabel.index = stationlabel.iloc[:,1]
             stationlabel.iloc[:,1] = stationlabel.iloc[:,0]
             
        #match station labels with lwcc labels
        label_lwcc = list(data.iloc[:,0])    
        label_overlap = list(set(label_station).intersection(label_lwcc))
        idx_label_station = [label_station.index(lb) for lb in label_overlap]
        idx_label_lwcc = [label_lwcc.index(lb) for lb in label_overlap]
        data = data.iloc[idx_label_lwcc,:]
        stationlabel = stationlabel.iloc[idx_label_station,:]

    
        data_salcorr = data.iloc[:,1:].values - stationlabel.iloc[:,-1].values[:,
                                    np.newaxis] * psi_sal.values 
    #    data_salcorr = pd.concat([stationlabel.iloc[:,1],pd.DataFrame(
    #            data_salcorr,index=stationlabel.index)], axis=1, sort=True)
        data_salcorr = pd.DataFrame(data_salcorr,index=stationlabel.index)
    #    data_salcorr.dropna(inplace=True)
        data_salcorr.sort_index(inplace=True)
        
        #smooth
        if smooth:
            figname = figname.replace('.png','_smooth.png')
            data_salcorr_smooth = savgol_filter(data_salcorr.interpolate(
                    limit_direction ='both'), 5, 1)
            data_salcorr = pd.DataFrame(data=data_salcorr_smooth, index=
                                        data_salcorr.index, columns=data_salcorr.columns)

        #null-point correction
        pos_NIR = np.where((wl>700) & (wl<720))[0]
        pos_null = np.where((wl>=690) & (wl<=710))[0]
        data_sal_nullcorr = data_salcorr
        for i in range(len(data_salcorr)):
            spec = data_salcorr.iloc[i,:]
            if any(spec.iloc[pos_NIR] < 0):
                data_sal_nullcorr.iloc[i,:] = spec - np.nanmedian(spec.iloc[pos_null].to_list())
            if any(spec.iloc[pos_NIR] > 0.05):
            #When the spectrum is significantly positive in the NIR (here set 
            #threshold as 0.05 m^-1), it is not considered as a valid measurement. 
                data_sal_nullcorr.iloc[i,:] = np.nan
        
        
        outname = figname.replace('.png','_null_final.txt')
        self.saveData(outname, data_sal_nullcorr, line)

    
        #plot        
        if plot:
            self.plotSpectrum(wl, data_sal_nullcorr.iloc[:,0:].transpose(), 
                              figname.replace('.png','_null_final.png'),
                              '$Salinity\ corrected\ a_{CDOM}$ [$m^{-1}$]')
    
            #Create target directory for individual spectral plots
            dir_fig = figname.replace('lwcc_merged_median_','').replace('.png','_null')
            self.createDir(dir_fig, overwrite=True)

            #plot ag spectra
            for j in range(len(data_sal_nullcorr)):           
                wl_for_plot = wl.iloc[wlpos['pos400']:wlpos['pos720']+1]
                data_for_plot = data_sal_nullcorr.iloc[j,wlpos['pos400']:wlpos['pos720']+1].to_numpy().transpose()
                namestr = os.path.join(dir_fig, data_sal_nullcorr.index[j] + '.png')
                self.plotSpectrum(wl_for_plot, data_for_plot,namestr, '$a_g(\lambda)\ [m^{-1}]$')
        
        print(f'Salinity effect corrected for {filename}!')
        
    #--------------------------------------------------------------------------
    def upload_pangaea(self, filename):
        '''
        This function restructs LWCC data suitable for uploading to Pangaea.
        Input:
        filename - str. Path of CDOM absorption or standard deviation data file (tab delimited).
        Output:
        LWCC data suitable to be uploaded to Pangaea saved as:
        "*_pangaea.txt", with "*" standing for the base name of input variable "filename".
        '''
        
        foo = filename.split('/')
        print(f'Importing {foo[-1]}...')
        
        wl, line, wlpos = self.getWL(filename)
        data = pd.read_csv(filename, header=None, comment='%', sep='\t')
        info = pd.read_csv(self.concurrent_info, comment='%', sep='\t')
        
        pos = pd.Series(np.nan, index=data.index)
        for i,lb in enumerate(data.iloc[:,0]):
            tmp = np.where(info.iloc[:,0]==lb)[0]
            if len(tmp) > 0:
                pos.iloc[i] = tmp[0]
                
        data_pangaea = info.iloc[pos,1:]
        data_pangaea['Date/Time'] = pd.to_datetime(data_pangaea['Date'] + ' ' + 
                                  data_pangaea['Time'])
        data_pangaea = data_pangaea.drop(['Date','Time'], axis=1)

        wl_for_data = wl.iloc[wlpos['pos350']:wlpos['pos730']+1]
        data_col = ['wl' + str(wl_for_data.iloc[i]) for i in range(len(wl_for_data))]
        new_data = pd.DataFrame(data=data.values[:,1+wlpos['pos350']:wlpos['pos730']+2], 
                                columns=data_col, index=data_pangaea.index)
        data_pangaea[data_col] = new_data
        
        #sort data by 'Date/Time'
        data_pangaea = data_pangaea.sort_values(by='Date/Time')
        
        #sort CTD data by Depth
        label_data_uni = np.unique(data_pangaea['Sample_ID'].tolist(),
                                   return_index=True, return_inverse=True, 
                                   return_counts=True)
        pos_ctd = np.where(label_data_uni[3] > 1)[0]
        for i in pos_ctd:
            pos = np.where(data_pangaea['Sample_ID']==label_data_uni[0][i])[0]
            tmp = data_pangaea.iloc[pos,:].sort_values(
                    by=['Depth_in_meter'])
            tmp.index = data_pangaea.iloc[pos,:].index
            data_pangaea.iloc[pos,:] = tmp
        
        new_col = ['Event', 'Date/Time', 'Sample_ID', 'Latitude', 'Longitude', 
                   'Depth_in_meter']
        data_pangaea = data_pangaea[new_col + [col for col in data_pangaea
                                                       if col not in new_col]]
        
        data_pangaea.to_csv(filename.replace('.txt','_pangaea.txt'), index=False, 
                           header=True, sep='\t', encoding='utf-8')
    
        print (f'Data from {foo[-1]} ready for Pangaea!')
    
    #--------------------------------------------------------------------------        
    def fit_exponential(self, filename, plot=True):
        '''
        This function fits CDOM absorption data with an exponential decay 
        function following Dall’Olmo et al. (2017), i.e. fitting region 420-490 nm,
        reference wavelength 440 nm.
        Input:
        1) filename - str. Path of CDOM absorption data file (tab delimited).
        2) plot - bool, default True. Plot fitted absorption spectra.
        Output:
        1) Exponentially fitted CDOM absorption data saved as 
            "*final_pangaea_expfit.txt". 
        2) Plots of the fitted data saved in the directory "*expfit".
        '''

        foo = filename.split('/')
        
        if 'pangaea' in filename:
            data = pd.read_csv(filename, sep='\t')
            pos = [col.find('wl') for col in data.columns].index(0)
            wl = pd.Series([float(str(wl).replace('wl','')) for wl in data.columns[pos:]])
            dir_fig = os.path.basename(filename).replace('lwcc_merged_median_', '').replace(
        '.txt', '') + '_expfit'
        else:
            data = pd.read_csv(filename, sep='\t', comment='%', header=None)
            wl, line, wlpos = self.getWL(filename)
            pos = 1
            dir_fig = os.path.basename(filename).replace('lwcc_merged_median_', '').replace(
        '.txt', '') + '_expfit' 
            
        wl_ref, wl_fit1, wl_fit2 = 440, 430, 490
        pos_fit1 = np.where(wl >= wl_fit1)[0][0]
        pos_fit2 = np.where(wl >= wl_fit2)[0][0]
        pos_null = np.where((wl>=690) & (wl<=710))[0]
        
        data_for_fit = data.iloc[:,pos+pos_fit1:pos+pos_fit2+1]
        data_for_fit = data_for_fit.dropna()
        wl_for_fit = wl.iloc[pos_fit1:pos_fit2+1]
        pos_ref = np.where(wl_for_fit >= wl_ref)[0][0]
        data_ref = data_for_fit.iloc[:,pos_ref]
        idx = data.index.intersection(data_for_fit.index)
        
        def exponential_func(x, slope, offset):
            #a_ref, wvl = x
            return x[0] * np.exp(slope * (wl_ref-x[1:]) ) + offset
        
        popt = [None] * len(data_for_fit)
        pcov = [None] * len(data_for_fit)
        
        data_esti = np.zeros(shape=(len(data_for_fit),len(wl)))
        
        #Create target directory for individual spectral plots         
        self.createDir(dir_fig, overwrite=True)
        
        for i in range(len(data_for_fit)):
            xdata = np.append(data_ref.values[i], wl_for_fit.values)
            popt[i], pcov[i] = curve_fit(exponential_func, xdata, data_for_fit.values[i,:])
            data_esti[i,:] = exponenial_func(np.append(data_ref.values[i],wl.values), *popt[i]) 
            data_esti[i,:] = data_esti[i,:] - np.nanmedian(data_esti[i,pos_null])      
            
            if plot:
                namestr = os.path.join(dir_fig, str(data.index[i]) + '.png')
                fig = plt.figure(figsize=(8.5, 6.5))
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.plot(wl, data_esti[i,:], 'r-',ls='--', label="Exponential Fit")
                plt.plot(wl, data.loc[idx,:].iloc[i,pos:], label = 'Salinity and null-corrected')
                plt.legend()
                fig.savefig(namestr, dpi=200)
                plt.close(fig)
                
        
        data_fit = pd.DataFrame(np.hstack([data.values[idx,:pos], data_esti]), 
                                columns=data.columns)
        if 'pangaea' in filename:
            self.saveData(filename.replace('.txt','_expfit.txt'), data_fit, 
                          index=False, header=True)
        else:
            self.saveData(filename.replace('.txt','_expfit.txt'), data_fit, 
                          headline=line, index=False)
    
        print (f'{foo[-1]} exponentially fitted!')

    #--------------------------------------------------------------------------        
    def abs_salcorr_pangaea(self):
        '''
        This function 
        1) extracts temperature and salinity data from the TSG data 
        downloaded from Pangaea, named as ”*phys_oce.tab“ and "*surf_oce.tab" 
        for CTD stations and underway data, respectively, for use (Note: If you 
        provide your own TS data, make sure they are named as 
        "*phys_oce_extracted.txt" and "*surf_oce_extracted.txt", respectively); 
        2) corrects salinity effects; 
        3) and restructs LWCC data suitable for uploading to Pangaea.
        '''
        
        if len(self.tsgST_pangaea)>0:
            filename_tsg = self.tsgST_pangaea.replace('.tab','_extracted.txt')
            if not os.path.isfile(filename_tsg):
                extract_pangaea_tempsal(self.tsgST_pangaea)
            filename = 'lwcc_merged_median_ag_station_uncorrected.txt'
            self.abs_salcorr0(filename, filename_tsg, self.labelST, 
                     smooth=True, plot=True)
            if self.TFpangaea:
                files = glob.glob('lwcc_merged*ag_station*.txt')
                for filename in files:
                    if ('pangaea' not in filename) & ('expfit' not in filename):
                        self.upload_pangaea(filename)

            
        if len(self.tsgUW_pangaea)>0:
            filename_tsg = self.tsgUW_pangaea.replace('.tab','_extracted.txt')
            if not os.path.isfile(filename_tsg):
                extract_pangaea_tempsal(self.tsgUW_pangaea, station=False)
            filename = 'lwcc_merged_median_ag_underway_uncorrected.txt'
            self.abs_salcorr0(filename, filename_tsg, self.labelUW, 
                     smooth=True, station=False, plot=True)
            if self.TFpangaea:
                files = glob.glob('lwcc_merged*ag_underway*.txt')
                for filename in files:
                    if ('pangaea' not in filename) & ('expfit' not in filename):
                        self.upload_pangaea(filename)

    
#-----------------------------------------------------------------------------
# Main Program. 
#-----------------------------------------------------------------------------
if __name__ == '__main__': 

    #Set working directory
    wd = '/work/ollie/yliu/Data/cruises/qftlwcc_postproc/lwcc'
    os.chdir(wd)

    process = lwccPostProc()
    
    #Step1:--------------------------------------------------------------------
    #process raw data.
    filenames_raw = sorted(glob.glob(os.path.join(wd,'*.raw')))
    for filename in filenames_raw:
        process.processRaw(filename)

    #Step2:--------------------------------------------------------------------
    #Plot CDOM absorption data from each sample file by file.  
    fnames_l0a = sorted(glob.glob(os.path.join(wd,'*.l0a')))
    for filename in fnames_l0a:
        process.plot_l0a(filename)

    #Step3:--------------------------------------------------------------------
    #1)remove spectra
    #Prepare your own 'rmSpectraIndex.txt' file (tab delimited) based on the 
    #plots generated from Step2 in the working directory!
    SpectraIndex = pd.read_csv('rmSpectraIndex.txt', comment='#', sep='\t')
    for i, filename in enumerate(SpectraIndex.iloc[:,0]):
        row_index = str(SpectraIndex.iloc[i,1])
        row_index = [float(idx) for idx in row_index.split(',')]
        process.rmabs_median(filename, row_index)
    
    #2)merge absorption and standard deviation data files and get all data lables.
    process.get_labels() 
    process.merge_abs()
    process.merge_abs(keyword='sd')
    
    #Step4:--------------------------------------------------------------------
    #1)Correct salinity effect and prepare data suitable for uploading to Pangaea
    # (whether the latter is needed is specified by the parameter "TFpangaea" in 
    #the configuration file "config_PostProc.txt").
    #Prepare your own salinity data for CTD stations and underway data, respectively.
    #If they are downloaded from Pangaea, they should be named as "*phys_oce.tab"
    #and "*surf_oce.tab", respectively; if not, rename the files as 
    #"*phys_oce_extracted.txt" and "*surf_oce_extracted.txt", respectively.
    #If data are to be uploaded to Pangaea, prepare your own "labels_latlon_datetime_lwcc.txt". 
    process.abs_salcorr_pangaea()
    
    #2)Fit CDOM absorption data with an exponential decay function following 
    #Dall’Olmo et al. (2017), i.e. fitting region 420-490 nm, reference 
    #wavelength 440 nm.
    files = glob.glob('*final*.txt')
    for filename in files:
        if 'expfit' not in filename:
            process.fit_exponential(filename, plot=True)
