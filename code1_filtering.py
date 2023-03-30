# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 09:57:04 2020

filtering background noise from sound files

@author: Noori Choi
"""
# Library imports
import os
import numpy as np
import itertools
import noisereduce as nr
import librosa
import librosa.display
from scipy.signal import butter, lfilter
from scipy.io import wavfile
from kneed import KneeLocator
from matplotlib import pyplot as plt

# Working directory
wd = "WORKING DIR"
# Listof soundfiles
sounds = os.listdir(wd)
# Directory where to save filtered file, must be different with wd
filt_loc = "DIR FOR SAVING FILTERED AUDIO FILES"
log = os.listdir(filt_loc)

plot = "recording plot ID"
date = "recording date"

# variables for bandpass filter
lowcut = 100
highcut = 5000
buffer = 0.1

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

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def bandpass_filter(buffer):
    return butter_bandpass_filter(buffer, lowcut, highcut, fs, order=5)

def find_silence(audio, fs, alpha, min_silence_sec=30, dt=0.01):
    """
    find the 'slience' in the audio file to detect sound event
    - modified from "pydub silence.py" to use adaptive threshold by sigma clipping
    - https://github.com/jiaaro/pydub/tree/master/pydub
    
    Sigma clipping methods:
        Amp_threshold = median(amplitude) + alpha*std(amplitude)
    
    parameters
    ----------
    audio : np.ndarray [shape=(# frames,) or (# channels, # frames)] 
        audio file to be segmented
    fs : int
        sampling rate
    alpha : float 
        optimal alpha for sigma clipping
    min_silence_sec: float (sec)
        minimum length of silence segment, default = 30
    dt: float (sec)
        minimum time interval for seek step, default = 0.01
        
    return
    ------
    silence_ranges : np.ndarray [shape = (2, # silence range)]
        begin and end of detected silence range (in sample number)
    """
    seg_len = len(audio)
    # Amplitude threshold
    ab_audio = np.absolute(audio)
    # set the seek_step and minimum noise length
    seek_step = int(fs * dt) 
    min_silence_len = int(fs*min_silence_sec) 
    
    # if audio is shorter than min_noise_len, return empty list
    if seg_len < min_silence_len:
        return []

    # setting silence threshold by sigma clipping
    silence_thresh = np.median(ab_audio) + (alpha*np.std(ab_audio))
    
    # find silence and add start and end indicies to the to_cut list
    silence_starts = []
    # check successive chunk of sound for silence
    # try a chunk at every "seek step"
    last_slice_start = seg_len - min_silence_len
    slice_starts = range(0, last_slice_start + 1, seek_step)

    # guarantee last_slice_start is included in the range
    # to make sure the last portion of the audio is searched
    if last_slice_start % seek_step:
        slice_starts = itertools.chain(slice_starts, [last_slice_start])

    for i in slice_starts:
        audio_slice = audio[i:i + min_silence_len]
        if np.max(audio_slice) <= silence_thresh:
            silence_starts.append(i)

    # short circuit when there is no silence
    if not silence_starts:
        return []

    # combine the silence we detected into ranges
    silence_ranges = []

    prev_i = silence_starts.pop(0)
    current_range_start = prev_i

    for silence_start_i in silence_starts:
        continuous = (silence_start_i == prev_i + seek_step)

        # sometimes two small blips are enough for one particular slice to be
        # non-silent, despite the silence all running together. Just combine
        # the two overlapping silent ranges.
        silence_has_gap = silence_start_i > (prev_i + min_silence_len)

        if not continuous and silence_has_gap:
            silence_ranges.append([current_range_start,
                                   prev_i + min_silence_len])
            current_range_start = silence_start_i
        prev_i = silence_start_i

    silence_ranges.append([current_range_start,
                           prev_i + min_silence_len])
    
    return silence_ranges

sounds = [i for i in sounds if f"{i[:-4]}_f.wav"  not in log]
for s in sounds:
    print(f"Processing {s}")
    #audio, fs = librosa.load(wd + '/' + s, sr=None, mono=True)
    fs, audio = wavfile.read(wd + '/' + s)
    #audio = highpass_filter(audio, 100, fs)
    
    # 1-1. Background noise filtering
    ## finding alpha for sigma clipping
    alpha = find_alpha(audio, fs, lowest=0, highest=5, num=20, knee_plot=False)
    ## extract silence range
    silence = find_silence(audio, fs, alpha, min_silence_sec=0.1, dt=0.01)
    if len(silence) == 0:
        ## if no silence detected, use the first 0.1 sec as the ref_silence
        print("no silence")
        silence_begin = 0
        silence_end = fs//10
        ref_silence = audio[silence_begin:silence_end] 
    else:
        ## find the longest silence
        ref_silence = silence[np.argmax([sil[1]-sil[0] for sil in silence])]
        silence_begin = ref_silence[0]
        silence_end = ref_silence[1]
        ref_silence = audio[silence_begin:silence_end]
    
    ## filtering raw audio file
    audio = nr.reduce_noise(y=audio, sr=fs, y_noise=ref_silence, 
                            n_fft=512, stationary=True)
    # 1-2. High-pass filter
    audio = np.apply_along_axis(bandpass_filter, 0, audio).astype('int16')
    wavfile.write(filt_loc + '/' + s[:-4] + '_f.wav', fs, audio)
    print("raw audio filtering is completed")

