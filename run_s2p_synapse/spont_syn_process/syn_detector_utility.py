#%%writefile MUST be the first line of the notebook to write it to a .py file

"""This utility function houses all the functions that we will use for correcting baseline fluorescence (BaselineRemoval /
 ZhangFit functions to remove photobleaching - mainly due to its use in whole-cell calcium imaging steps I do.
 PEAK_SETTINGS are not used for synapses and are dealt with in the single_spine functions listed below)"""

# This will write a python file to be used in other files (e.g. from spont_syn_detector_utility import * == aka from this file
# import ALL functions)
import pandas as pd
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

from scipy.signal import find_peaks, peak_prominences
from BaselineRemoval import BaselineRemoval
# from oasis.functions import deconvolve



"""
Settings han; determine empirically later you can largely ignore these unless you are using the bottom detector functions (oasis / hatem)
"""
# PEAK_SETTINGS = {
#     "distance": 20,
#     "height": 7,
# } LEFT OVER FROM MARTI ESTABLISHING ARBITRARY CUTOFFS


def single_spine_peak_plotting(input_f, input_fneu):
    corrected = input_f - (0.7*input_fneu)
        # a standard subtraction of fluorescence background immediately surrounding each ROI
    corrected_sample_trace = BaselineRemoval(corrected)
    corrected_trace = corrected_sample_trace.ZhangFit(repitition = 1000, lambda_ = 500)
        # an adaptively weighted iterated modified polynomial fit that ignores peaks and corrects baseline to 0
    mini_peaks, _ = find_peaks(corrected, height = 2*(abs(np.median(corrected_trace)) + abs(corrected_trace.min())), distance = 40)
        # scipy find_peaks function
        #then plot the traces you generate
    plt.plot(corrected)
    plt.plot(mini_peaks, corrected[mini_peaks], "x")
    plt.plot(np.full_like(corrected, threshold), "--",color = "grey")
    plt.plot(np.full_like(corrected, np.median(corrected)), "--", color = 'r')
    plt.show()
    
""" https://suite2p.readthedocs.io/en/latest/outputs.html explains the code to load a video's F / Fneu files (and others)
    With the F / Fneu files; we can iterate this function to see the outputs using the following:
    for f, fneu in zip(F, Fneu):
        single_spine_peak_plotting(f, fneu)"""


def single_synapse_baseline_correction_and_peak_return(input_f, input_fneu, return_peaks = True):
    
    corrected_trace = input_f - (0.7*input_fneu)
    corrected_trace = BaselineRemoval(corrected_trace)
    corrected_trace = corrected_trace.ZhangFit(repitition = 100)
    peaks, _ = find_peaks(corrected_trace, height = 2*(abs(np.median(corrected_trace)) + abs(corrected_trace.min())), distance = 40)
    amplitudes = corrected_trace[peaks] - np.median(corrected_trace)
    
    if return_peaks == True:
        return peaks
    else:
        return amplitudes

""" please go from this file to the suite2p_utility notebook to see how this will be used within Suite2p's output """ 

def fluorescence_decay():
    """Here we need to determine the number of frames after a peak"""    
    
    
 #The following are  lines of code written by Marti, who I generally trust to write cohesive and concise code
# for coding; If you find my functions above to be inadequate you may impliment these instead
    
    

def detect_spikes_by_mod_z(input_trace, **signal_kwargs):
    median = np.median(input_trace)
    deviation_from_med = np.array(input_trace) - median
    mad = np.median(np.abs(deviation_from_med))
    mod_zscore = deviation_from_med/(1.4826*mad)
    return signal.find_peaks(mod_zscore, **signal_kwargs)[0]


def plot_spikes(raw_trace, detector_func, detector_trace=None, **detector_kwargs):
    if detector_trace is None:
        detector_input_trace = raw_trace.copy()
    else:
        detector_input_trace = detector_trace.copy()
    spikes = detector_func(detector_input_trace, **detector_kwargs)
    plt.plot(range(len(raw_trace)), raw_trace, color="blue")
    for spk in spikes:
        plt.axvline(spk, color="red")
    plt.show()
    
    
    
    
# rolling_min / remove_bleaching are basic polynomial-based baseline corrections; this will not remove noise either
# these can be compared to the more complicated ZhangFit iteratively-weighted approach (airPLS method)

def rolling_min(input_series, window_size):
    r = input_series.rolling(window_size, min_periods=1)
    m = r.min()
    return m


def remove_bleaching(input_trace):
    min_trace = rolling_min(pd.Series(input_trace), window_size=int(len(input_trace)/10))
    fit_coefficients = np.polyfit(range(len(min_trace)), min_trace, 2)
    fit = np.poly1d(fit_coefficients)
    return input_trace - fit(range(len(input_trace)))
