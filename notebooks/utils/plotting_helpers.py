# -*- coding: utf-8 -*-
"""
Helper functions used to simplify code for visualization of results. 
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# GRABBING FUNCTIONS FOR CREATING PLOT DFS
def grab_polymer_name(full_filepath, common_filepath):
    '''Grabs the polymer name from a file path.
    
    Parameters:
    -----------
    full_filepath: string
        path to the datasheet for that polymer 

    common_filepath: string 
        portion of the filepath that is shared between all polymer inputs, excluding their custom names
    
    Returns:
    -------
    polymer_name: string
        the custom portion of the filepath with the polymer name and any other custom info
    '''

    _, polymer_name = full_filepath.split(common_filepath)
    polymer_name = polymer_name[:-5]  # remove the .xlsx

    return polymer_name

def grab_peak_subset_data(df, ppi):
    '''Subsets a df to only the desired proton peak index.

    Parameters:
    ----------
    df: Pandas.DataFrame
        replicate all data for one polymer ("all" flag signals raw data is from all protons, not just binding ones)
    ppi: int
        the desired proton peak index to subset to

    Returns:
    -------
    subset_df: Pandas.DataFrame
        copy of the dataframe after subset to the desired peak
    '''
    subset_df = df.loc[df['proton_peak_index'] == ppi].copy()

    return subset_df

def grab_disco_effect_data(subset_df):
    '''Grabs basic statistical summary of each disco effect(t) 
    datapoint for a given peak in a given polymer.

    Parameters:
    -----------
    subset_df: Pandas.Dataframe
        subset_df comes is the desired peak subset of the replicate all raw data for one polymer
    
    Returns:
    --------
    grouped_df: pandas.DataFrame
        snapshot of all mean statistical parameters for each sat time and selected peak

    mean_disco: pd.Series
        mean disco effect value for the given peak at all sat times
    
    std_disco: pd.Series
        std dev disco effect value for the given peak at all sat times
    
    n: int
        number of replicates associated with statistical summary 
    '''

    grouped_df = subset_df.groupby(by=['sat_time', 'proton_peak_index']).mean()
    mean_disco = subset_df.groupby(by=['sat_time', 'proton_peak_index']).mean()['corr_%_attenuation']
    std_disco = subset_df.groupby(by=['sat_time', 'proton_peak_index']).std()['corr_%_attenuation']
    n = subset_df['replicate'].max()

    return grouped_df, mean_disco, std_disco, n

def assemble_peak_buildup_df(df, ppi):
    '''Assembles a directly indexable dataframe of individual proton Disco Effect(t) data, 
    that can be used to plot an individual peak buildup curve for the desired proton.

    Parameters:
    -----------
    df: pandas.Dataframe
        replicate all data file from one polymer, contains data for all peaks

    ppi: int
        proton peak index of the desired polymer proton for the buildup curve
    
    Returns:
    --------
    plot_df: pandas.DataFrame
        the assembled df containing data required to plot the buildup curve of the selected peak

    '''
    # grab data for desired peak
    subset_df = grab_peak_subset_data(df, ppi)

    # manipulate raw data and reassemble into a flat single column index plot df
    grouped_df, _, std_disco, n = grab_disco_effect_data(subset_df)
    plot_df = grouped_df.copy()
    plot_df = plot_df.rename(columns={"corr_%_attenuation": "corr_%_attenuation_mean"})
    plot_df['corr_%_attenuation_std'] = std_disco
    plot_df['sample_size'] = n
    plot_df = plot_df.drop(columns='replicate') # replicate data column doesn't make sense after grouping by mean
    plot_df = plot_df.reset_index()

    return plot_df

# ANNOTATION FUNCTIONS
def annotate_sig_buildup_points(ax, significance, sat_time, disco_effect, dx, dy, color):
    '''Adds markers to the buildup curve plot provided for by ax where points are flagged as
    statistically significant in the passed "significance" series. 
    
    Parameters:
    ----------
    ax: matplotlib plot axis
        indicates the plot to be annotated

    significance: pd.Series or array-like
        iterable the length of the plot domain containing boolean flags for domain indices that should
        or should not be annotated as significant
    
    sat_time: np.array
        plot domain values

    disco_effect: np.array
        plot range values
    
    dx: float
        amount to shift the annotation marker away from datapoint on x axis

    dy: float:
        amount to shift the annotation marker away from datapoint on y axis

    color: string 
        color-code or other matplotlib compatible signifier of marker colour

    Returns:
    -------
    None, performs annotation action on the plot ax

    '''

    sig_annotation_markers = significance.map({True: "*", False: " "}).reset_index(drop=True)

    for ix, marker in enumerate(sig_annotation_markers):
        ax.annotate(marker, xy=(
            sat_time[ix]+dx, disco_effect[ix]+dy), c=color)

    return

