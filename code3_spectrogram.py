# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 14:43:16 2022

1. load files
2. pad the signals (if it's needed)
3. extracting log spectrogram from signals
4. normalize spectrogram
5. save the normalized spectrogram

PreprocessingPipeline

@author: Noori
"""
import librosa
import numpy as np
import os
import pickle
from scipy.signal import butter, filtfilt

class Loader:
    """loader is responsible for loading on an audio file"""
    
    def __init__(self, sample_rate, frame_size, hop_length, frequency_threshold,\
                 silence_threshold, mono):
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_length = hop_length
        self.frequency_threshold = frequency_threshold
        self.silence_threshold = silence_threshold
        self.mono = mono
    
    def load(self, file_path):
        signal = librosa.load(file_path, 
                              sr = self.sample_rate,
                              mono = self.mono)[0]
        return signal
    
    def lowpass_filter(self, signal, sample_rate, frequency_threshold):
        cutoff = frequency_threshold/(0.5*sample_rate)
        b, a = butter(9, cutoff, btype='low', analog=False)
        signal = filtfilt(b, a, signal)
        return signal
    
    def trim_silence(self, signal, sample_rate, frame_size, hop_length, silence_threshold):
        """
        trim silence at the beginning and end of the audio file
        
        parameters
        ----------
        signal : np.ndarray [shape=(# frames,) or (# channels, # frames)] 
            raw audio file to be trimmed
        sample_rate : float
            sampling rate of the input audio file
        frame_size : int
            frame length for rms calculation
        hop_length : int
            hop length for rms calculation
        silence_threshold : float (sec)
            minimum amplitude threshold for determining silence
            proportion to the maximum rms value (min-max normalized)
        
        returns
        -------
        signal : np.ndarray [shape=(# frames,) or (# channels, # frames)] 
            trimmed audio file
        """
        # calculate and normalize rms values of audio file
        signal_rms = librosa.feature.rms(y=signal, frame_length=frame_size, 
                                        hop_length=hop_length)[0]
        signal_rms = MinMaxNormalizer.normalize(signal_rms)
        
        # find silence at the beginning
        trim_start = 0 # ms
        while signal_rms[trim_start:trim_start+1] < silence_threshold and trim_start < len(signal_rms):
            trim_start += 1
        
        # find silence at the end
        trim_end = -1 # ms
        while signal_rms[trim_end-1:trim_end] < silence_threshold and -1*trim_end < (len(signal_rms)+1):
            trim_end += -1
        
        # convert rms frame into frame in audio file
        trim_start = trim_start*hop_length
        trim_end = trim_end*hop_length
        
        # trim silence from raw audio file
        signal = signal[trim_start:trim_end]
        return signal
    
class Padder:
    """Padder is responsible to apply padding to an array"""
    
    def __init__(self, mode="constant"):
        self.mode = mode
        
    def left_pad(self, array, num_missing_items):
        padded_array = np.pad(array,
                              (num_missing_items, 0),
                              mode = self.mode)
        return padded_array
    
    def right_pad(self, array, num_missing_items):
        padded_array = np.pad(array,
                              (0, num_missing_items),
                              mode = self.mode)
        return padded_array
    
class LogSpectrogramExtractor:
    """LogSpectrogramExtractor extracts log spectrogram (dB) from a time-series signal"""
    
    def __init__(self, frame_size, hop_length):
        self.frame_size = frame_size
        self.hop_length = hop_length
        
    def extract(self, signal):
        stft = librosa.stft(signal,
                            n_fft = self.frame_size,
                            hop_length = self.hop_length)[:-1]
        spectrogram = np.abs(stft)
        log_spectrogram = librosa.amplitude_to_db(spectrogram)
        return log_spectrogram

class MinMaxNormalizer:
    """MinMaxNormalizer applies Min-Max normalization to an array"""
    
    def __init__(self, min_val, max_val):
        self.min = min_val
        self.max = max_val
        
    def normalize(self, array):
        norm_array = (array - array.min()) / (array.max() - array.min())
        norm_array = norm_array * (self.max - self.min) + self.min
        return norm_array
    
    def denormalize(self, norm_array, original_min, original_max):
        array = (norm_array - self.min) / (self.max - self.min)
        array = array * (original_max - original_min) + original_min
        return array

class Saver:
    """saver is responsible to save features, and min max values"""
    
    def __init__(self, feature_save_dir, min_max_values_save_dir, plot):
        self.feature_save_dir = feature_save_dir
        self.min_max_values_save_dir = min_max_values_save_dir
        self.plot = plot
        
    def save_feature(self, feature, file_path):
        save_path = self._generate_save_path(file_path)
        np.save(save_path, feature)
        return save_path
        
    def save_min_max_values(self, min_max_values):
        save_path = os.path.join(self.min_max_values_save_dir,
                                 f"{self.plot}_min_max_values.pkl")
        self._save(min_max_values, save_path)
    
    @staticmethod
    def _save(data, save_path):
        with open(save_path, "wb") as f:
            pickle.dump(data, f)
    
    def _generate_save_path(self, file_path):
        file_name = os.path.split(file_path)[1]
        save_path = os.path.join(self.feature_save_dir, file_name[:-4] + ".npy")
        return save_path      

class PreprocessingPipeline:
    """PreprocessingPipeline processes audio files in a directory, 
    applying the following steps to each file"""
    
    def __init__(self):
        self.padder = None
        self.extractor = None
        self.normalizer = None
        self.saver = None
        self.min_max_values = {}
        self._loader = None
        
    @property
    def loader(self):
        return self._loader
    
    @loader.setter
    def loader(self, loader):
        self._loader = loader
        
    def process(self, audio_files_dir):
        for root, _, files in os.walk(audio_files_dir):
            for file in files:
                file_path = os.path.join(root, file)
                self._process_file(file_path)
                print(f"Processed file {file_path}")
        self.saver.save_min_max_values(self.min_max_values)
    
    def _process_file(self, file_path):
        signal = self.loader.load(file_path)
        
        feature = self.extractor.extract(signal)
        norm_feature =  self.normalizer.normalize(feature)
        save_path = self.saver.save_feature(norm_feature, file_path)
        self._store_min_max_value(save_path, feature.min(), feature.max())
        
    def _store_min_max_value(self, save_path, min_val, max_val):
        self.min_max_values[save_path] = {
            "min": min_val,
            "max": max_val
            }
        
if __name__ == "__main__":
    FRAME_SIZE = 512
    HOP_LENGTH = 256
    LOW_PASS = 8000 # in Hz
    SILENCE_THRESHOLD = 0.05 # in seconds
    SAMPLE_RATE = None
    MONO = True
    
    plot_date_list = ["PLOT1_DATE1", "PLOT2_DATE2"]
    
    for PLOT in plot_date_list:
        SPECTROGRAM_SAVE_DIR = "DIR FOR SAVING SPECTROGRAMS"
        MIN_MAX_VALUES_SAVE_DIR = "DIR FOR SAVIMG MIN_MAX_VALUES"
        FILES_DIR = "WORKING DIR"
        
        # instantiate objects
        loader = Loader(SAMPLE_RATE, FRAME_SIZE, HOP_LENGTH, LOW_PASS, SILENCE_THRESHOLD, MONO)
        padder = Padder()
        log_spectrogram_extractor = LogSpectrogramExtractor(FRAME_SIZE, HOP_LENGTH)
        min_max_normalizer = MinMaxNormalizer(0, 1)
        saver = Saver(SPECTROGRAM_SAVE_DIR, MIN_MAX_VALUES_SAVE_DIR, PLOT)
        
        preprocessing_pipeline = PreprocessingPipeline()
        preprocessing_pipeline.loader = loader
        preprocessing_pipeline.padder = padder
        preprocessing_pipeline.extractor = log_spectrogram_extractor
        preprocessing_pipeline.normalizer = min_max_normalizer
        preprocessing_pipeline.saver = saver
        
        preprocessing_pipeline.process(FILES_DIR)
    
    
        



