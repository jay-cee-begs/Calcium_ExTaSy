
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

frame_rate = 10
EXPERIMENT_DURATION = 180
#total number of seconds that the experiment lasts (tends to be flexible so we will have to figure out how to set this)
FRAME_INTERVAL = 0.1
#frame_interval is calculated as 1 / frame_rate
#we are not doing binned experiments so this is not necessary; again, will comment out later, but for now too scared of bugs
BIN_WIDTH = 10
#FILTER_NEURONS applies strictly to 'iscell.npy' file; in most instances, we will use all ROIs anyway, but keep true for clarity
FILTER_NEURONS = True



def list_all_files_of_type(input_path, filetype):
    return [file for file in os.listdir(input_path) if file.endswith(filetype)]

def string_to_list_translator(input_string, strip_before_split="[ ]", split_on=" "):
    split_string = input_string.strip(strip_before_split).split(split_on)
    return list(filter(None, split_string))

def spike_list_translator(input_string):
    """This funciton is nested in the next. It is designed to convert the time stamp of each event into a time
        during the experiment (e.g. frame 2 = 1.1 seconds into the recording)"""
    string_list = string_to_list_translator(input_string)
    return np.array(string_list).astype(int) * FRAME_INTERVAL

def amplitude_list_translator(input_string):
    amp_string_list = string_to_list_translator(input_string)
    amp_string_list = np.array(amp_string_list).astype(float)
    return np.around(amp_string_list)

def decay_list_translator(input_string):
    decay_string_list = string_to_list_translator(input_string)
    decay_string_list = np.array(decay_string_list).astype(int)
    return np.array(decay_string_list).astype(int)
    #need to find a way to calculate the difference between decay time and decay points; decide which one we are going to save in the pipeline

def spike_df_iterator(input_path, return_name=True):
    for csv_file in list_all_files_of_type(input_path, "csv"):
        csv_path = os.path.join(input_path, csv_file)
        csv_df = pd.read_csv(csv_path, converters={"PeakTimes":spike_list_translator , "Amplitudes":amplitude_list_translator})
        yield csv_df, csv_file if return_name else csv_df

        
   #Again, binned stats are not necessary for synapses, but regardless, we can leave this here for now
def calculate_binned_stats(input_df):
    local_df = input_df.copy()

    bins = np.arange(0, EXPERIMENT_DURATION + BIN_WIDTH, BIN_WIDTH)
    population_spikes, _ = np.histogram(np.hstack(local_df["PeakTimes"].values), bins=bins)
    population_frequency = population_spikes / BIN_WIDTH

    bin_stats = pd.DataFrame.from_dict({
        "Bin_Limits": [(bins[bin_index], bins[bin_index + 1]) for bin_index in range(len(bins) - 1)],
        "Spikes": population_spikes,
        "Frequency": population_frequency})
        
    return bin_stats


def calculate_cell_freq(input_df):
    output_df = input_df.copy()
    output_df["SpikesCount"] = output_df["PeakTimes"].str.len()
    output_df["SpikesFreq"] = output_df["SpikesCount"] / (input_df["Total Frames"] * frame_rate) #divide by total # of frames NOT framerate
    return output_df

def calculate_cell_isi(input_df): #isi == interspike interval
    output_df = input_df.copy()
    output_df["SpikesDiff"] = output_df["PeakTimes"].apply(lambda x: list(pd.Series(x).diff().dropna()))
    output_df["DiffAvg"] = output_df["SpikesDiff"].apply(lambda x: pd.Series(x).mean())
    output_df["DiffMedian"] = output_df["SpikesDiff"].apply(lambda x: pd.Series(x).median())
    output_df["DiffCV"] = output_df["SpikesDiff"].apply(lambda x: pd.Series(x).std()) / output_df["DiffAvg"] * 100
    return output_df

#below I will need to accurately figure out how to integrate this in, I should meet with Marti Ritter next week to do so

def calculate_spike_amplitudes(input_df):
    output_df = input_df.copy()
    output_df["AvgAmplitude"] = output_df["Amplitudes"].apply(lambda x: pd.Series(x).mean())
    output_df["SpkAmpMedian"] = output_df["Amplitudes"].apply(lambda x: pd.Series(x).median())
    output_df["SpkAmpCV"] = output_df["Amplitudes"].apply(lambda x: pd.Series(x).std()) / output_df["AvgAmplitude"] * 100
    return output_df


def calculate_cell_stats(input_df, calculate_freq=True, calculate_isi=True, calculate_amplitudes=True):
    output_df = input_df.copy()
    if calculate_freq:
        output_df = calculate_cell_freq(output_df)
    if calculate_isi:
        output_df = calculate_cell_isi(output_df)
    if calculate_spike_amplitudes:
        output_df = calculate_spike_amplitudes(output_df)
    return output_df



def process_spike_csvs_to_pkl(input_path, output_path, overwrite=False):
    """This will convert .csv files into pickle files which behave like dataframes; but are faster and preserve CPU RAM"""
    for spike_df, file_name in spike_df_iterator(input_path):
        processed_path = os.path.join(output_path, 
                                      f"{os.path.splitext(file_name)[0]}"
                                      f"Dur{int(EXPERIMENT_DURATION)}s"
                                      f"Int{int(FRAME_INTERVAL*1000)}ms"
                                      f"Bin{int(BIN_WIDTH*1000)}ms"
                                        + ("_filtered" if FILTER_NEURONS else "") +
                                      ".pkl")

        if os.path.exists(processed_path) and not overwrite:
            print(f"Processed file {processed_path} already exists!")
            continue
            
        if FILTER_NEURONS:
            spike_df = spike_df[spike_df["IsUsed"]]
            
        processed_dict = {
            "cell_stats": calculate_cell_stats(spike_df),
            "binned_stats": calculate_binned_stats(spike_df)}
        

        pd.to_pickle(processed_dict, processed_path)
