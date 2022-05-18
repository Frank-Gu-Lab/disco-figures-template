# -*- coding: utf-8 -*-
"""
Custom scripts to plot DISCO profiles to a passed plot axis. 
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import numpy as np
from utils.plotting_helpers import annotate_sig_buildup_points
from utils.wrangle_data import flatten_multicolumns, calculate_abs_buildup_params

# set default style params
plt.style.use(['science', 'discolib'])
plt.rcParams.update({'font.family': 'sans-serif'})
plt.rcParams.update({'font.size': 6})


# for use in more specific plotting
custom_colors = ['#377eb8', '#984ea3', '#ff7f00', '#e41a1c', '#f781bf',
    '#ffff33', '#4daf4a', '#a65628', '#999999']

# must install LaTex before you can use Science Plots
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'


def add_fingerprint_toax(df, ax, **kwargs):
    '''Adds DISCO AF0 fingerprints of all binding protons for one polymer to the plot axis passed.
        
    Parameters:
    -----------
    df: pandas.DataFrame
        "replicate" raw data for one polymer (raw file containing all technical AF0 replicates for only the polymers binding protons)
    ax: matplotlib axis
        axis of the plot to contain the DISCO fingerprint
    
    **kwargs: dict {"custom_palette": [list, of, colourcodes]} (optional, if desired)
        custom plotting palette containing a unique color code for each proton 

    Returns:
    --------
    None, adds DISCO AF0 fingerprint to plot axis passed

    Example:
    -------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt

    >>> # read the per-replicate data table, containing information from only the binding protons in HPC
    >>> hpc_df_replicates = pd.read_excel("../data/raw/stats_analysis_output_replicate_HPC_370k_20uM.xlsx", index_col=[0], header=[0]).reset_index()
    >>> fig, ax = plt.subplots(figsize=(5, 5))
    >>> add_fingerprint_toax(df=hpc_df_replicates, ax=ax)
    >>> plt.show()

    The binding AF0 fingerprint for HPC will appear on the passed plot axis.
    '''

    ppm = np.round(df.copy().groupby(by='proton_peak_index')['ppm'].mean(), 2)
    ppm_ix = np.unique(df['proton_peak_index'].values)

    # map ppm to proton peak index incase multi ppms per proton peak index
    ppm_mapper = dict(zip(ppm_ix, ppm))
    df['ppm'] = df['proton_peak_index'].map(ppm_mapper)

    # remove duplicates from df
    df_plot = df[['ppm', 'AFo']].drop_duplicates()

    # take absolute value of AFo
    df_plot['AFo'] = df_plot['AFo'].abs()

    #create plot, use custom palette if there is one
    try:
        custom_palette = kwargs.pop("custom_palette")
        sns.barplot(data=df_plot, x='ppm', y='AFo', ax=ax,
                    edgecolor='k',errcolor='black',palette=custom_palette)
        sns.stripplot(data=df_plot, x='ppm', y='AFo', ax=ax, edgecolor='k',
                      linewidth=0.5, jitter=False, palette=custom_palette)

    except KeyError: # executes if there is no custom palette
        pass
        sns.barplot(data=df_plot, x='ppm', y='AFo', ax=ax, edgecolor='k', errcolor='black')
        sns.stripplot(data=df_plot, x='ppm', y='AFo', ax=ax,
                      edgecolor='k', linewidth=0.5, jitter=False)

    # reformat x axis
    ax.invert_xaxis()  # invert to match NMR spectrum

    return

def add_buildup_toax(df, ax):
    '''Adds the Absolute DISCO Effect(t) buildup curves for all protons in 
    a polymer to an existing plot axis.
    
    Parameters:
    -----------
    df: pandas.DataFrame
        mean raw data file for one polymer, containing mean DISCO Effects(t) 
        (corr % attuenation) for /only/ the binding protons 

    ax: matplotlib axis
        plot axis to contain the buildup curves
    
    Returns:
    --------
    None, adds buildup curves to axis

    Example:
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pandas as pd
    >>> # read the mean data table, containing information from only the binding protons in HPC
    >>> hpc_df = pd.read_excel("../data/raw/stats_analysis_output_mean_HPC_370k_20uM.xlsx", index_col=[0, 1, 2, 3], header=[0, 1]).reset_index()
    >>> fig, ax = plt.subplots(figsize = (5,5))
    >>> add_buildup_toax(df = hpc_df, ax = ax)
    >>> plt.show()
    
    The buildup curves for all binding peaks in HPC will appear on the plot axis.
    '''

    if type(df.columns) == pd.MultiIndex:
        df = flatten_multicolumns(df)  # make data indexable if it is not already

    ppm_labels = np.round(df['ppm'], 2)
    df_plot = df.copy()
    groups = df_plot.groupby([ppm_labels])

    # plot DISCO effect build up curve, absolute values
    for ppm, group in groups:
        sat_time, disco_effect, y1, y2 = calculate_abs_buildup_params(group)
        ax.plot(sat_time, disco_effect, markeredgecolor='k', markeredgewidth=0.35,
                marker='o', linestyle='', ms=5, label="%.2f" % ppm) 
        ax.fill_between(sat_time, y1, y2, alpha=0.25)
        ax.legend(loc='lower right', title="Peak ($\delta$, ppm)")
        ax.axhline(y=0.0, color="0.8", linestyle='dashed')

    return

def add_overlaid_buildup_toax_customlabels(df_list, ax, **kwargs):
    '''Adds buildup curve(s) to a plot. Allows control of individual proton buildup curves when designing a buildup curve plot, 
    for use cases where only specific subsets of buildup curves are desired to be plotted.
    
    Parameters:
    -----------
    df_list: list
        list of one "mean" raw dataframe for every polymer whose data is to be used in the plot
    
    ax: matplotlib axis
        location for the plot
    
    **kwargs: dict 
        any desired custom plotting parameters

    Notes:
    ------
    * To add new custom properties, simply add another kwargs.pop statement below, and pass a new kwarg to the dict
    
    Example:
    -------
    >>> import matplotlib.pyplot as plt 
    >>> import pandas as pd 
    >>> from utils.plotting_helpers import assemble_peak_buildup_df 
    >>> from utils.wrangle_data import generate_disco_effect_mean_diff_df, generate_subset_sattime_df

    >>> # read data
    >>> low_hpc_replicate_all = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_80k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
    >>> high_hpc_replicate_all = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_370k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
    
    >>> # grab the desired proton to plot by their proton peak index
    >>> ppi_1_low = assemble_peak_buildup_df(low_hpc_replicate_all, 1)
    >>> ppi_1_high = assemble_peak_buildup_df(high_hpc_replicate_all, 1) 

    >>> # compute change significance between proton buildup curves
    >>> hpc_effect_size_df = generate_disco_effect_mean_diff_df(low_hpc_replicate_all, high_hpc_replicate_all)
    >>> hpc_subset_sattime_df = generate_subset_sattime_df(hpc_effect_size_df, 0.25)
    >>> kwargs = {"labels": ["80", "370"],
    >>>       "dx": 0.003,
    >>>       "dy": 0.010,
    >>>       "change_significance": hpc_effect_size_df,
    >>>       "annot_color": "#000000",
    >>>       "custom_colors": ['#b3cde3', '#377eb8']}

    >>> # construct and display peak overlay plot w/ annotated change significance
    >>> fig, ax = plt.subplots(figsize = (5,2))
    >>> df_list = [ppi_1_low, ppi_1_high]
    >>> add_overlaid_buildup_toax_customlabels(df_list, ax, **kwargs)
    >>> plt.show()
    '''

    # extract custom properties
    custom_labels = kwargs.pop("labels")
    dx = kwargs.pop("dx")
    dy = kwargs.pop("dy")
    change_sig_df = kwargs.pop("change_significance")
    buildup_colors = kwargs.pop("custom_colors")
    annot_color = kwargs.pop("annot_color")

    # plot overlaid buildups using the correct custom properties
    color_count = 0
    for ix, df in enumerate(df_list):
        plot_label = custom_labels[ix]

        if type(df.columns) == pd.MultiIndex:
            df = flatten_multicolumns(df)  # make data indexable if it is not already

        ppm_labels = np.round(df['ppm'], 2)
        df_plot = df.copy()
        groups = df_plot.groupby([ppm_labels])

        # plot DISCO effect build up curve, absolute values
        for _, group in groups:
            sat_time, disco_effect, y1, y2 = calculate_abs_buildup_params(group)

            full_plot_label = f"{plot_label}"
            ax.plot(sat_time, disco_effect, markeredgecolor='k', markeredgewidth=0.35, color=buildup_colors[color_count],
                    marker='o', linestyle='', ms=5, label=full_plot_label)
            ax.fill_between(sat_time, y1, y2, color=buildup_colors[color_count],
                            alpha=0.25)
            ax.axhline(y=0.0, color="0.8", linestyle='dashed')
        color_count += 1

    # annotate significance of change in disco effect (NOT disco adhesion interaction significance)
    key = group['proton_peak_index'].unique()[0]
    change_sig_subset = change_sig_df.loc[change_sig_df['proton_peak_index'] == key]

    # annotate change sig points
    significance = change_sig_subset['changed_significantly']
    annotate_sig_buildup_points(ax, significance, sat_time, disco_effect, dx, dy, color=annot_color)

    return

def add_difference_plot(df, ax, dy, **kwargs):
    '''Add a proton-wise difference profile plot to a plot axis (Vertical Orientation). 
    
    Parameters:
    -----------
    df: Pandas.DataFrame
        "subset_sattime_df" for a polymer comparison, contains change effect data computed only at
        the desired sat time t for this plot
    
    ax: matplotlib axis
        plot destination

    dy: float
        y distance for significance marker away from datapoint

    **kwargs: dict {"custom_colors":[list, with, one, custom, color, per, proton]} (optional)

    Returns:
    --------
    None, adds difference profile to a plot
    
    Theory Notes:
    -------------
    * Proton-wise difference profile plot shows the change effect size and significance of any changes in mean DISCO effect(t) between all
    the protons from two polymer test groups. In this work, test groups were low mW vs high mW versions of the same polymer.

    * Selected sat time for the difference plot is determined when the subset_sattime_df is computed. Lower sat times are more
    representative as there is the least influence from rebinding.

    Example:
    --------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from utils.plotting import add_difference_plot
    >>> from utils.wrangle_data import generate_disco_effect_mean_diff_df, generate_subset_sattime_df

    >>> # we use "replicate all" data as we are considering both binding and non-binding protons
    >>> low_HPC = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_80k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
    >>> high_HPC = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_370k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
        
    >>> # first compute HPC's difference profile from the raw data
    >>> hpc_effect_size_df = generate_disco_effect_mean_diff_df(low_HPC, high_HPC) 

    >>> # subset HPC's difference profile to the desired sat time Disco Effect(t = 0.25)
    >>> hpc_subset_sattime_df = generate_subset_sattime_df(hpc_effect_size_df, 0.25) 

    >>> kwargs = {"custom_colors": ['#377eb8', '#984ea3', '#ff7f00', '#e41a1c', '#f781bf', '#ffff33', '#4daf4a', '#a65628', '#999999']}

    >>> fig, ax = plt.subplots(figsize = (3,5))
    >>> add_difference_plot(df=hpc_subset_sattime_df, dy= 0.010,  ax=ax, **kwargs)
    >>> plt.show()
    '''
    custom_colors = kwargs.pop("custom_colors")

    plot_range = range(1, (df.shape[0])+1)

    # zero line
    ax.axvline(x=0.0, color="0.8", linestyle='dashed')

    # error bars
    ax.hlines(y=plot_range, xmin=df['effect_sem_lower'],
              xmax=df['effect_sem_upper'], color='black', linewidth=2, zorder=1)
 
    # data
    ax.scatter(df['effect_size'], plot_range, s=(40,), color=custom_colors[:df.shape[0]],
               alpha=1, label='Effect Size', marker='o', linewidths=0.35, edgecolors='k', zorder = 2)
    ax.set_yticks(plot_range, np.round(df['ppm'].values, 2))
    
    # annotate significance
    df['annotation'] = df['changed_significantly'].map({True: "*", False: ""})

    for ix, value in enumerate(list(plot_range)):
        x = df['effect_size'].iloc[ix]
        y = value + dy
        marker = df['annotation'].iloc[ix]
        ax.annotate(marker, (x, y), c='#000000')

    return

def add_difference_plot_transposed(df, ax, dy, **kwargs):
    '''Add a proton-wise difference profile plot to a plot axis (Horizontal Orientation). 
    
    Parameters:
    -----------
    df: Pandas.DataFrame
        "subset_sattime_df" for a polymer comparison, contains change effect data computed only at
        the desired sat time t for this plot
    
    ax: matplotlib axis
        plot destination

    dy: float
        y distance for significance marker away from datapoint

    **kwargs: dict {"custom_colors":[list, with, one, custom, color, per, proton]} (optional)

    Returns:
    --------
    None, adds difference profile to a plot
    
    Theory Notes:
    -------------
    * Proton-wise difference profile plot shows the change effect size and significance of any changes in mean DISCO effect(t) between all
    the protons from two polymer test groups. In this work, test groups were low mW vs high mW versions of the same polymer.

    * Selected sat time for the difference plot is determined when the subset_sattime_df is computed.

    Example:
    --------
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from utils.plotting import add_difference_plot_transposed
    >>> from utils.wrangle_data import generate_disco_effect_mean_diff_df, generate_subset_sattime_df

    >>> # we use "replicate all" data table as we are considering both binding and non-binding protons
    >>> low_HPC = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_80k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
    >>> high_HPC = pd.read_excel("../data/raw/stats_analysis_output_replicate_all_HPC_370k_20uM.xlsx", index_col=[0], header=[0]).reset_index(drop=True)
        
    >>> # first compute HPC's difference profile from the raw data
    >>> hpc_effect_size_df = generate_disco_effect_mean_diff_df(low_HPC, high_HPC) 

    >>> # subset HPC's difference profile to the desired sat time Disco Effect(t = 0.25)
    >>> hpc_subset_sattime_df = generate_subset_sattime_df(hpc_effect_size_df, 0.25) 

    >>> kwargs = {"custom_colors": ['#377eb8', '#984ea3', '#ff7f00', '#e41a1c', '#f781bf', '#ffff33', '#4daf4a', '#a65628', '#999999']}

    >>> fig, ax = plt.subplots(figsize=(3, 5))
    >>> add_difference_plot_transposed(df=hpc_subset_sattime_df, dy=0.010, ax=ax, **kwargs)
    >>> plt.show() 
    '''

    custom_colors = kwargs.pop("custom_colors")

    plot_domain = range(1, df.shape[0]+1)

    # zero line
    ax.axhline(y=0.0, color="0.8", linestyle='dashed')
    
    # error bars
    ax.vlines(x=plot_domain, ymin=df['effect_sem_lower'],
              ymax=df['effect_sem_upper'], color='black', linewidth=2, zorder = 1)

    # data
    ax.scatter(plot_domain, df['effect_size'], s = (40,), color=custom_colors[:df.shape[0]],
               alpha=1, label='effect size', marker = 'o', linewidths = 0.35, edgecolors = 'k', zorder = 2)    
    ax.set_xticks(plot_domain, np.round(df['ppm'].values,2))

    # annotate significance
    df['annotation'] = df['changed_significantly'].map({True: "*", False: ""})
    
    for ix, value in enumerate(list(plot_domain)):
            y = df['effect_size'].iloc[ix] + dy
            x = value + 0.05 
            marker = df['annotation'].iloc[ix]
            ax.annotate(marker, (x,y), c = '#000000')

    return

