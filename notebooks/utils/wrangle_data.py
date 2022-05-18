# -*- coding: utf-8 -*-
"""
Data processing scripts involved in preprocessing, transforming, and statistically testing data for plotting.
"""

import numpy as np
import pandas as pd 
import scipy.stats as stats
from scipy.stats import t
import pingouin as pg

# DATA PREPROCESSING FUNCTIONS
def flatten_multicolumns(mean_df):
    '''Takes in a mean df and flattens the multi-tier column
    index into a directly indexable index.'''

    # clean up multi index for both
    colnames = mean_df.columns.get_level_values(0).values
    mean_df = mean_df.droplevel(1, axis=1)
    colnames[4] = "corr_%_attenuation_mean"
    colnames[5] = "corr_%_attenuation_std"
    mean_df.columns = colnames
    
    return mean_df

# CALCULATE DISCO EFFECT PARAMS
def calculate_abs_buildup_params(df):
    ''' Calculates Absolute buildup curve parameters from a plot-level group
    of data.

    Parameters:
    ----------
    df: pd.DataFrame
        mean corr % attenuation and associated data for a proton-level buildup curve
        after being flattened by flatten_multicolumns, sourced originally from mean raw data files

    Returns:
    --------
    sat_time: array-like
        saturation times for the buildup curve

    disco_effect: array-like
        disco effect for the buildup curve

    y1: array-like
        standard error lower bound 

    y2: array-like
        standard error upper bound
    '''

    # plot DISCO effect build up curve, absolute values
    sat_time = df['sat_time'].values
    disco_effect = abs(df['corr_%_attenuation_mean'].values)
    std = abs(df['corr_%_attenuation_std'].values)
    n = df['sample_size'].values
    std_err = std/np.sqrt(n)

    y1 = np.subtract(disco_effect, std_err)
    y2 = np.add(disco_effect, std_err)

    return sat_time, disco_effect, y1, y2

def calculate_buildup_params(df):
    ''' Calculates buildup curve parameters (not absolute) from a plot-level group
    of data.

    Parameters:
    ----------
    df: pd.DataFrame
        corr % attenuation and associated data for a proton-level buildup curve
        after being flattened by flatten_multicolumns

    Returns:
    --------
    sat_time: array-like
        saturation times for the buildup curve

    disco_effect: array-like
        disco effect for the buildup curve

    y1: array-like
        standard error lower bound 

    y2: array-like
        standard error upper bound
    '''

    # plot DISCO effect build up curve, absolute values
    sat_time = df['sat_time'].values
    disco_effect = df['corr_%_attenuation_mean'].values
    std = df['corr_%_attenuation_std'].values
    n = df['sample_size'].values
    std_err = std/np.sqrt(n)

    y1 = np.subtract(disco_effect, std_err)
    y2 = np.add(disco_effect, std_err)

    return sat_time, disco_effect, y1, y2

# CHANGE PROFILE DISCO EFFECT STAT TESTING, AND RESULTING DATA WRANGLING FUNCTIONS
def shapiro_wilk(effect):
    '''
    H0: the data is normally distributed
    H1: the data is not normally distributed
    
    Therefore, if p > 0.05 we fail to reject the null hypothesis'''

    stat, p = stats.shapiro(effect)

    return p

def bartlett(effect1, effect2):
    '''After testing that p(shapiro wilk) > 0.05 for both effect1 and effect2, use Bartlett's test to
    examine the equal variance assumption.

    H0: the data is of equal variance 
    H1: the data is not of equal variance
    
    Therefore, if p > 0.05 we fail to reject the null hypothesis
    '''

    stat, p = stats.bartlett(effect1, effect2)

    return p

def test_change_significance(group_1, group_2, alt_hyp="two-sided"):
    '''Defaults to two sided test, can change to 'greater'  
    for seeing if effect is larger in group 2 than group 1.
    or 'less' for vice versa.

    GROUP 1 = CONTROL (LOW)
    GROUP 2 = TREATMENT (HIGH)

    Reports back the significance of any changes, and their effect size for interpretation.
    '''
    df_list = []

    for first, second, in zip(group_1, group_2):

        # assign data and indexes of each group included in statistical testing
        first_ix = first[0]
        first_df = first[1]
        second_ix = second[0]
        second_df = second[1]

        if first_ix == second_ix:  # validates that the same sat time and proton peak are being compared btw two groups
            ppm = first_df['ppm'].values[0]
            ppi = first_df['proton_peak_index'].values[0]

            disco_1 = abs(first_df['corr_%_attenuation'])
            disco_2 = abs(second_df['corr_%_attenuation'])

            # test for normality
            norm_p1 = shapiro_wilk(disco_1)  # if fail to reject null, they are treated as normal
            norm_p2 = shapiro_wilk(disco_2)

            # test for equal variance
            p_eq_var = bartlett(disco_2, disco_1)

            # if norm dist and equal variance do parametric test
            if np.logical_and(np.logical_and(norm_p2, norm_p1), p_eq_var) > 0.05:

                # calc and flag significance results
                stat, p_result = stats.ttest_ind(disco_2, disco_1, equal_var=True, alternative=alt_hyp)

                # calc 95 confidence intervals of datapoints
                # python code ref: https://stats.stackexchange.com/questions/475289/confidence-interval-for-2-sample-t-test-with-scipy
                disco_1_mean = np.mean(disco_1)
                v1 = np.var(disco_1, ddof = 1)
                n1 = len(disco_1)

                disco_2_mean = np.mean(disco_2)
                v2 = np.var(disco_2, ddof = 1)
                n2 = len(disco_2)
                
                delta_disco = disco_2_mean - disco_1_mean
                pooled_se = np.sqrt(v1 / n1 + v2 / n2)
                dof = (v1 / n1 + v2 / n2)**2 / (v1**2 / (n1**2 * (n1 - 1)) + v2**2 / (n2**2 * (n2 - 1)))
                
                # upper and lower bounds 95 CI
                delta_CI_lower = delta_disco - t.ppf(0.95, dof)*pooled_se 
                delta_CI_upper = delta_disco + t.ppf(0.95, dof)*pooled_se

                if float(p_result) <= 0.05:
                    sig_bool = True
                    print(f"Sig Point is: {first_ix[0]}, {ppm}, p = {p_result}, n = {len(disco_1)}")

                else:
                    sig_bool = False

                # calc effect size
                hedges_g = pg.compute_effsize(disco_2, disco_1, paired=False, eftype="hedges")
                
                # calc effect size sem for error bars
                # see pingouin docs for formula https://pingouin-stats.org/generated/pingouin.compute_esci.html
                effect_se = np.sqrt(((n1+n2) / (n1*n2)) + ((hedges_g**2) / (2*(n1+n2)))) 

                # write results to output
                current_dict = {"sat_time": first_ix[0],
                                "proton_peak_index": ppi,
                                "ppm": ppm,
                                "changed_significantly": sig_bool,
                                "mean_difference":delta_disco,
                                "delta_95_ci_lower": delta_CI_lower,
                                "delta_95_ci_upper": delta_CI_upper,
                                "delta_sem_lower": delta_disco - pooled_se,
                                "delta_sem_upper": delta_disco + pooled_se,
                                "effect_size": hedges_g,
                                "effect_sem_lower":hedges_g - effect_se,
                                "effect_sem_upper":hedges_g + effect_se}

                df_list.append(current_dict)

            else:
                print("Observations do not pass normality and equal variance test!")
                print("Shapiro P's: ", norm_p1, norm_p2)
                print("Equal Variance P:", p_eq_var)

    results_df = pd.DataFrame(df_list)

    return results_df

def generate_subset_sattime_df(effect_size_df, sat_time):
    '''Subset the disco observations to only those of the desired sat_time, to generate
    the desired subset data for the difference profile plot.'''

    subset_sattime_df = effect_size_df.loc[effect_size_df['sat_time'] == sat_time].copy().drop(columns="sat_time")

    return subset_sattime_df

def generate_disco_effect_mean_diff_df(replicate_df_low, replicate_df_high):
    '''Take raw dfs from input file and generate dfs
    contining information about the change in disco effect (mean difference) between them.
    Effect size (hedge's g), effect significance with 95% CI.
    
    Notes:
    ------
    Must be dfs from two identical polymers with the same peaks. I.e.
    low mW and high mW. Note that the convention for ordering matters.

    In all calculation cases, order is [high df - low df], as is low df is the control
    and high df is the treatment. High df is the df with the property being 
    increased, i.e. increased mW or % hydrolysis.
    '''

    # take absolute values of everything
    grouped_low = replicate_df_low.groupby(
        by=["sat_time", "proton_peak_index"])
    grouped_high = replicate_df_high.groupby(
        by=["sat_time", "proton_peak_index"])

    # perform t test per peak per sat time to see if sig change w increased property
    effect_size_df = test_change_significance(group_1=grouped_low, group_2=grouped_high)

    return effect_size_df

