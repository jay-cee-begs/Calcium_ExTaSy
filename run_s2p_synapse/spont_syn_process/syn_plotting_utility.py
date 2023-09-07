#I probably will rename these files so I can tell them apart when I have a lot of tabs open...

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import os

"""
below is an example structure for a dictionary for all the experiment files 
suite2one might implement this code as the following in a processing_pipeline:
    experiment_structure = {
    "control": {
        "replicate_1": [f"control_A_well_0{j}Dur180sInt100msBin500ms_filtered.pkl" for j in range(1, 5)],
        "replicate_2": [f"control_A_well_0{j}Dur180sInt100msBin500ms_filtered.pkl" for j in range(5, 9)],
    },
    "treatment": {
        "replicate_1": [f"treatment_B_well_0{j}Dur180sInt100msBin500ms_filtered.pkl" for j in range(1, 5)],
        "replicate_2": [f"treatment_B_well_0{j}Dur180sInt100msBin500ms_filtered.pkl" for j in range(5, 9)],
    }
}
experiment_structure = {

"APV_treatment": [f"DIV1{i}_cs{j}_r{k}_APV*.pkl" for i in range(5,9), j in range(1,3), k in range (1,3)],
"PDBu_treatment": [f"DIV1{i}_cs{j}_r{k}_PDBu*.pkl" for i in range(5,9), j in range(1,3), k in range (1,3)],
"CNQX_treatment": [f"DIV1{i}_cs{j}_r{k}_CNQX*.pkl" for i in range(5,9), j in range(1,3), k in range (1,3)],
"baseline": [f"DIV1{i}_cs{j}_r{k}_baseline*.pkl" for i in range(5,9), j in range(1,3), k in range (1,3)]

maybe consider a dictionary of dictionaries?? although I am not sure the benefit of this immediately
}
"""
_experiment_structure_example = {
    "control": {
        "dataset1": ["file1", "file2"]
    },
    "APV": {
        "dataset2": ["file3", "file4"]
    },
    "PDBu": {
        "dataset3": ["file1", "file2"]
    },
    "CNQX": {
        "dataset4": ["file3", "file4"]
    },
}

_available_tests = {
    "mann-whitney-u": stats.mannwhitneyu,
    "wilcoxon": stats.wilcoxon,
    "paired_t": stats.ttest_rel,
}
#If you cannot tell this is where I have spent the least time so far in coding...

def get_significance_text(series1, series2, test="mann-whitney-u", bonferroni_correction=1, show_ns=False, 
                          cutoff_dict={"*":0.05, "**":0.01, "***":0.001}, return_string="{text}\n~{pvalue:.4f}"):
    statistic, pvalue = _available_tests[test](series1, series2)
    levels, cutoffs = np.vstack(list(cutoff_dict.items())).T
    levels = np.insert(levels, 0, "n.s." if show_ns else "")
    text = levels[(pvalue < cutoffs.astype(float)).sum()]
    return return_string.format(pvalue=pvalue, text=text)

def add_significance_bar_to_axis(ax, series1, series2, center_x, line_width):
    significance_text = get_significance_text(series1, series2, show_ns=True)
    
    original_limits = ax.get_ylim()
    
    ax.errorbar(center_x, original_limits[1], xerr=line_width/2, color="k", capsize=4)
    ax.text(center_x, original_limits[1], significance_text, ha="center", va="bottom")
    
    extended_limits = (original_limits[0], (original_limits[1] - original_limits[0]) * 1.2 + original_limits[0])
    ax.set_ylim(extended_limits)
    
    return ax

def aggregated_feature_plot(experiment_df, feature="SpikesFreq", agg_function="median", comparison_function="mean",
                            palette="Set1", significance_check=["control", "treatment"]):
    
    grouped_df = experiment_df.groupby(["dataset", "group", "file_name"]).agg(agg_function).reset_index(drop=False)

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot()
    color_palette = sns.color_palette(palette)

    sns.swarmplot(x="group", y=feature, hue="dataset", data=grouped_df, ax=ax, size=10, alpha=0.6, palette=palette)

    marker_width = 0.5
    tick_positions = {ax.get_xticklabels()[index].get_text(): ax.get_xticks()[index] for index in
                      range(len(ax.get_xticklabels()))}
    
    for group, group_df in grouped_df.groupby("group"):
        for dataset_index, (dataset, dataset_df) in enumerate(group_df.groupby("dataset")):
            feature_mean = dataset_df[feature].agg(comparison_function)
            ax.plot([tick_positions[group] - marker_width/2, tick_positions[group] + marker_width/2],
                    [feature_mean, feature_mean], "-", color=color_palette[dataset_index], lw=2)

    if significance_check:
        sub_checks = [significance_check] if not any(isinstance(element, list) for element in significance_check) else significance_check
        for sub_check in sub_checks:
            add_significance_bar_to_axis(ax, 
                                 grouped_df[grouped_df["group"] == sub_check[0]][feature], 
                                 grouped_df[grouped_df["group"] == sub_check[1]][feature],
                                (tick_positions[sub_check[0]] + tick_positions[sub_check[1]]) / 2,
                                 abs(tick_positions[sub_check[0]] - tick_positions[sub_check[1]]))

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, frameon=False)

    ax.set_ylim([0, ax.get_ylim()[1]])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig

def build_experiment_dfs(input_path, experiment_structure):
    experiment_cell_stats, experiment_binned_stats = pd.DataFrame(), pd.DataFrame()

    for group in experiment_structure.keys():
        for dataset in experiment_structure[group].keys():
            for file_name in experiment_structure[group][dataset]:
                file_dict = pd.read_pickle(os.path.join(input_path, file_name))
                cell_stats, binned_stats = file_dict["cell_stats"], file_dict["binned_stats"]
                
                for stats, experiment_stats in zip((cell_stats, binned_stats), 
                                                   (experiment_cell_stats, experiment_binned_stats)):
                    stats[["group", "dataset", "file_name"]] = group, dataset, file_name
                experiment_cell_stats = pd.concat((experiment_cell_stats, cell_stats))
                experiment_binned_stats = pd.concat((experiment_binned_stats, binned_stats))
    return experiment_cell_stats, experiment_binned_stats
        
