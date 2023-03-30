# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:51:01 2020

peakfind & bout grouping with GMM

@author: Noori Choi
"""
# Library imports
# Library imports
import os
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn import mixture
import librosa
import librosa.display
from scipy.io import wavfile

plot = "recording plot name"
date = "date of audio recording"
# Working directory
wd = "WORKING DIR"
# Listof soundfiles
sounds = os.listdir(wd)
# Directory where to save the csv file, must be different with wd
csv_loc = "DIR FOR SAVING DETECTED PEAKS DATAFRAME"
nopulse_loc = "DIR FOR SAVING THE LIST OF SILENT AUDIO FILES"
seg_loc = "DIR FOR SEGMENTED DETECTED AUDIO FILES"
loglist = os.listdir(csv_loc)

buffer = 0.1 # buffers for audio segmentation (in seconds)

def find_alpha(audio, fs, lowest=0, highest=10, num=20, 
               alpha_plot=False, knee_plot=False):
    """
    find the optimal alpha for sigma clipping methods by knee(elbow) method
    
    Sigma clipping methods:
        Amp_threshold = median(amplitude) + alpha*std(amplitude)
    
    parameters
    ----------
    audio : np.ndarray [shape=(# frames,) or (# channels, # frames)] 
        audio file to be segmented
    fs : int
        sampling rate
    lowest : float 
        lowest alpha, default = 0
    highest : float
        highest alpha, default = 10
    num : int
        steps for np.linspace(), default = 20
    alpha_plot : boolean
        visualize sigma clipping method, default = False
    knee_plot : boolean
        save knee(elbow) plot, default = False
    
    returns
    -------
    knee : float
        optimal alpha value found by knee(elbow) method
    """
    ab_audio = np.absolute(audio)
    # create the list of alpha values
    list_alpha = np.linspace(lowest, highest, num)
    
    n_sample = []
    for idx, alpha in enumerate(list_alpha):
        #set amp_threshold
        A_low = np.median(ab_audio) + (alpha*np.std(ab_audio))
        #count number of samples below amp_threshold
        sample = np.count_nonzero(audio > A_low)
        n_sample.append([alpha, sample])
        
        if alpha_plot == True:
            # make a figure
            fig, ax = plt.subplots(figsize=(14,10))
            librosa.display.waveshow(y=audio, sr=fs, x_axis='time', ax=ax)
            ax.set_title('waveform', fontweight="bold", size=20)
            ax.tick_params(axis='both', which='major', labelsize=20)
            ax.set_xlabel('Time (sec)', fontsize = 22)
            ax.axhline(A_low, ls='-', color='r')
            plt.show()
            #plt.savefig(csv_loc + f"/alpha_{s}_{idx}.png")
            plt.close()
           
    n_sample = np.array(n_sample) 
    
    # knee(elbow) method
    ## may need to check the knee plot for the first time 
    ## to set parameters (curve, direction)
    kn = KneeLocator(n_sample[:, 0], n_sample[:, 1], curve='convex', direction='decreasing')
    if knee_plot == True:
        kn.plot_knee()
        plt.show()
        plt.close()
        #plt.savefig(filt_loc + f'/kn_plot_{s}')
    
    if kn.knee:
        knee = kn.knee
    else:
        raise ValueError("cannot find kneed point in the given range of alphas")
    return knee

def peak_find(file_name, audio, alpha, dt=0.03):
    """
    find the time point of peaks above amplitude threshold in audio files
    
    Amplitude threshold by Sigma clipping methods:
        Amp_threshold = median(amplitude) + alpha*std(amplitude)
    
    parameters
    ----------
    file_name : string
        name of audio file
    audio : np.ndarray [shape=(# frames,) or (# channels, # frames)] 
        audio file to be segmented
    alpha : float
        alpha for sigma clipping method
    dt : float
        the minimum interval between two peaks (in seconds), default = 0.03 (s)
    
    returns
    -------
    result_df : dataframes
        dataframes with time and amplitude of detected peaks
    """
    # sigma clipping for peak detection threshold
    A_low = np.median(audio) + (alpha*np.std(audio))
    # detect peaks
    peaks, properties = find_peaks(audio, height = A_low, distance = dt*fs)
    time = np.ndarray.tolist(peaks/fs) # time of peaks in second
    amp = np.ndarray.tolist(properties["peak_heights"]) # amp of each peak
    peaks = np.column_stack((time, amp))
    # create dataframe
    result = np.hstack((peaks, np.array([[file_name]] * len(peaks))))
    result_df = pd.DataFrame({'Time': result[:, 0], 'Amp': result[:, 1], 'File': result[:,2]})
    return result_df

def feature_extraction(audio, fs, n_fft=2048, hop_length=512):
    """
    Calculate spectral features for excluding noise clips

    Parameters
    ----------
    audio : np.ndarray
        Audio time series, Multi-channel is supported.
    fs : Number > 0 [scalar]
        frame rate
    n_fft : int > 0, optional
        FFT window size. The default is 2048.
    hop_length : int > 0, optional
        Hop length for STFT. The default is 512.

    Returns
    -------
    feature_ls : list
        list of calculated values
    """
    flat = librosa.feature.spectral_flatness(audio, n_fft=n_fft, hop_length=hop_length)[0]
    bw = librosa.feature.spectral_bandwidth(audio, sr=fs, n_fft=n_fft, hop_length=hop_length)[0]
    rms = librosa.feature.rms(y=audio, frame_length=n_fft, hop_length=hop_length)[0]
    
    feature_ls = [max(flat), min(flat), np.median(flat), 
                  max(bw), min(bw), np.median(bw),
                  max(rms), min(rms), np.median(rms)]
    return feature_ls

for index, item in enumerate(loglist):
    loglist[index] = loglist[index][:-4]

loglist = []
no_pulse = []
too_noisy = []
for s in sounds:
    if s[:-4] in loglist:
        print(s + " is already processed")
        pass
    else:
        audio, fs = librosa.load(wd + '/' + s, sr=None)
        # Find peaks in each file
        ## get alpha for sigma clipping
        alpha = find_alpha(audio, lowest=0, highest=3, num=10)
        ## find peaks using determined alpha
        df = peak_find(s, audio, alpha, dt=0.03)
        if len(df) < 20:
            ## if there are only a few peaks, remove the file from further analysis
            no_pulse.append(s)
        elif len(df) > (len(audio)//(3*fs))//0.03:
            ## if there are too many peaks, remove the file from further analysis
            too_noisy.append(s) 
        else:
            # Group detected peaks into a signal bout
            ## calculate time interval
            df['Time'] = pd.to_numeric(df['Time'])
            df['dT'] = df['Time'].diff().fillna(0)
            df = df[1:]
            df['dT'] = np.log(df['dT'])
            
            df_index = df.reset_index(drop = True)
            df_dt = df[['dT']]
            X = df_dt.to_numpy().reshape(-1, 1)
            ## Gaussian Mixture Model
            high_thr = np.log(1) # set outlier max limit
            low_thr = np.log(0.03)
            interval = df[['dT']].loc[(df.dT > df.dT.quantile(.05)) & (df.dT < df.dT.quantile(.95))]
            ### find best_fit
            N = np.arange(3, 11)
            models = [None for i in range(len(N))]
            
            for i in range(len(N)):
                models[i] = mixture.GaussianMixture(N[i]).fit(interval)
            
            #### compute the AIC and the BIC
            best_fit = []
            for j in range(25):
                AIC = [m.aic(X) for m in models]
                BIC = [m.bic(X) for m in models]
                bout_gen = models[np.argmin(AIC)]
                best_fit.append(bout_gen)
            
            bout_gen = max(set(best_fit), key=best_fit.count)
            bout_gen.fit(X) #GMM model fit
            probs = bout_gen.predict_proba(X) #Soft clustering
            
            probs_df = pd.DataFrame(data=probs, columns=[str(i) for i in range(probs.shape[1])])
            probs_df['group'] = probs_df.idxmax(axis=1) #assign dT into the highest probability
            fin_df = pd.concat([df_index, probs_df], axis=1)
            
            bout_df = fin_df.groupby('group')['dT'].median()
            criteria_b = bout_df.idxmax()
            
            ## Bout grouping by the largest group of dT
            fin_df['bout'] = (fin_df['group'] == criteria_b).groupby(fin_df['File']).cumsum() + 1
            
            ## Pulse grouping by GMM
            fin_df = fin_df.groupby(['File', 'bout']).agg(
                begin = ('Time', min),
                end = ('Time', max)).reset_index()
            
            
            fin_df['duration'] = fin_df['end'] - fin_df['begin']
            fin_df = fin_df[fin_df['duration'] > buffer*2]
            
            # 3. Audio segmentation
            time_stamp = fin_df[['bout', 'begin', 'end']].to_numpy()
            
            feature_ls = []
            for n_bout in range(len(time_stamp)):
                bout_name = int(time_stamp[n_bout][0])
                bout_begin = int(max(time_stamp[n_bout][1]*fs - buffer*fs, 0))
                bout_end = int(min(time_stamp[n_bout][2]*fs + buffer*fs, len(audio)))
                ## extract spectral features (bandwidth & flatness)
                audio_seg = audio[bout_begin:bout_end]
                spec_features = feature_extraction(audio_seg, fs)
                feature_ls.append(spec_features)
                ## generate audio segment between adjacent silences
                wavfile.write(seg_loc + '/' + s[:-4] + f'_{bout_name}.wav', fs, audio_seg.astype('int16'))
                
            ## append spectral features to ind_df
            spec_columns = ['max_flat', 'min_flat', 'med_flat', 
                            'max_bw', 'min_bw', 'med_bw',
                            'max_rms', 'min_rms', 'med_rms']
            fin_df[spec_columns] = np.array(feature_ls)
            ## export csv file with time stamp
            fin_df.to_csv(csv_loc + f'/time_stamp_{s[:-4]}.csv', index=False, sep = ",")
            
        no_pulse_df = pd.DataFrame({'file':no_pulse})
        no_pulse_df.to_csv(nopulse_loc + f'/{date}_nopulse.csv', index=False, sep=",")
        
        noisy_df = pd.DataFrame({'file':too_noisy})
        noisy_df.to_csv(nopulse_loc + f'/{date}_noisy.csv', index=False, sep=",")
        