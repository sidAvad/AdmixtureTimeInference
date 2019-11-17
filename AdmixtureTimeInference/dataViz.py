import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ast
from numpy import random
import seaborn as sns

# Pre-processing functions

def consolidate_data(input_directory, pedsim_postfix, t, ntypes):
    """
    Data for a sigle admixture type is consolidated across all simulation batches
    Returns a dictionary of dataframes, with one entry per admixture type

    Args:
    input_directory - Input directory containing results
    t - Admixture time simulated ( used to locate input filepaths)
    nsims - numbers of simulation batches
    ntypes - number of admixture types

    Returns:
    dataframe - consolidated dataframe
    """
#{{{ 
    appended_data_types = []
    for admixtype in range(ntypes):

        # Put all the different simulation batches together
        input_file = input_directory + '/{}gens/admixture-time_{}gens_type{}_{}.res'.format(t, t, admixtype + 1,
                                                                                            pedsim_postfix)

        current_dataframe = pd.read_table(input_file)
        current_dataframe['admixtype'] = admixtype + 1

        appended_data_types.append(current_dataframe)
#}}}

    return (pd.concat(appended_data_types, ignore_index=True))


def merge_string_lists(row, cols):
    '''
    Evaluates and merges two lists as strings. Returns merged list. To be used with pandas df.apply
    '''
    # {{{

    # for i in range(len(cols)):
    #    returnlist = returnList + row[cols].tolist()  
    #    print(returnlist)

    returnlist = row[cols[0]] + row[cols[1]] + row[cols[2]]
    # }}}
    return (returnlist)

def get_mle(target, time, typ, m, finite_chromosome_correction=True):
    '''
    Takes in unphased diploid list of tract lengths ,admixture time, admixture type and m count which is the number of segements that end with a chromosome end rather than a breakpoint. Returns mle for admixture time. 
    '''
    #{{{ 
    if len(target) == 0:
        mle = None
    else:
        mle = (len(target) - m) / (sum(target) * 2) if finite_chromosome_correction else len(sample_list) / (sum(sample_list) * 2)
    #}}}

    return (mle)


# ----------Functions for final mle estimation based on simone gravel model---------------#

def get_mle_list_twopart(dataframe, target_hom1, target_hom2, time, typ, mcount_column,
                         finite_chromosome_correction=True):
    '''
    Takes in dataframe and returns mle estimates for target columns ( het1 and het2; should contain lists of unphased diploid tract lenghts) for each row(sample) given admixture time, admixture type and m counts which is a column containing number of segements that end with a chromosome end rather than a breakpoint.This is the version of the function that estimates mle for a subset of tract types
    '''
    # {{{

    subset_dataframe = dataframe[(dataframe.admixtype == typ) & (dataframe.t == time)]
    tract_lengths_list1 = list(subset_dataframe[target_hom1])
    tract_lengths_list2 = list(subset_dataframe[target_hom2])
    mle_list = []

    for indx, sample_list in enumerate(tract_lengths_list1):
        m = subset_dataframe[mcount_column].iloc[indx]

        mle_linearcomb1 = (len(sample_list1) - m) / (sum(sample_list1)) + 1 if finite_chromosome_correction else len(
            sample_list1) / (sum(sample_list)) + 1
        mle_linearcomb2 = (len(sample_list2) - m) / (sum(sample_list2)) + 1 if finite_chromosome_correction else len(
            sample_list2) / (sum(sample_list)) + 1

        mles = get_simultaneous_mle(mle_linearcomb1, mle_linearcomb2, p1, p2)

        mle_list.append(mles)

    # }}}
    return (mle_list)


def get_simultaneous_mle(rhs0, rhs1, p1, p2):
    # {{{
    if p1 == 0 or p1 == 1:
        t2 = 2
        t1 = 2
    else:
        t2 = (p1 - p1(rhs0 + rhs1) + rhs1 - 1) / (p1 - p1 * p2 + p2) + 1
        t1 = (rhs0 - (1 - p2) * t2) / (1 - p1)
    # }}}
    return (t1, t2)


def get_simult_mle_df(hom0_list, hom1_list, p1, p2, m_hom0, m_hom1, mode = "default"):

    # {{{

    """
    Function to get simult mle for use with df
    """ 

    try:
        rhs0 = (len(hom0_list) - m_hom0) / (sum(hom0_list))
        rhs1 = (len(hom1_list) - m_hom1) / (sum(hom1_list))
    except:
        print('homozygous tract lists empty; returning NA')
        return (None, None)

    if mode == "joint":
        return((rhs0 + rhs1 +2)/2)

    if p1 == 0 or p1 == 1:
        t2 = 2
        t1 = 2
    else:
        A = np.array([[p1, p2], [1 - p1, 1 - p2]])
        b = np.array([rhs1, rhs0])
        t = np.linalg.solve(A, b) + 1
    

    # }}}
    return t


def get_simult_mle_diagnose_df(hom0_list, hom1_list, p1, p2, m_hom0, m_hom1, mode = "default"):

    # {{{

    """
    Function to get simult mle for use with df
    """ 

    try:
        rhs0 = (len(hom0_list) - m_hom0) / (sum(hom0_list))
        rhs1 = (len(hom1_list) - m_hom1) / (sum(hom1_list))
    except:
        print('homozygous tract lists empty; returning NA')
        return (None, None)

    if mode == "joint":
        return((rhs0 + rhs1 +2)/2)

    if p1 == 0 or p1 == 1:
        t2 = 2
        t1 = 2
    else:
        A = np.array([[p1, p2], [1 - p1, 1 - p2]])
        b = np.array([rhs1, rhs0])
        t = np.linalg.solve(A, b) + 1
    
    # }}}

    return [sum(hom0_list),sum(hom1_list), rhs0, rhs1]

def compute_markovWF(p1,p2,t1,t2):
    return [p1*(t1-1)+p2*(t2-1), (1-p1)*(t1-1) + (1-p2)*(t2-1)]

# Plotting functions ( single )
def plot_histograms_side_types(dataframe, targets, admixture_types):
    '''
    Plots side-by-side histograms of target column/s (containing lists) of dataframe given list ofadmixture types to plot for
    '''
    # {{{
    tract_lengths_all_types = []
    for i, typ in enumerate(admixture_types):

        if isinstance(targets, str):
            tract_lengths = dataframe[dataframe.admixtype == typ][
                targets]  # Subset entries to admixture_type of interest
        else:
            tract_lengths = dataframe[dataframe.admixtype == typ][
                targets[i]]  # Subset entries to admixture_type of interest

        tract_lengths_list = list(tract_lengths)
        tract_lengths_all_sims = [x for y in tract_lengths_list for x in y]  # putting together all simulations
        tract_lengths_all_types.append(tract_lengths_all_sims)

    # plot list of lists of tractlenths ( one for each admixture type) as a side by side histogram
    plt.hist(tract_lengths_all_types, bins=20, density=True, label=admixture_types)
    plt.xlabel('tract lengths')
    plt.ylabel('Frequency')
    plt.title('tract length distributions')
    # }}}
    return ()


def plot_histograms_side(dataframe, target, dictionary):
    '''
    Plots side-by-side histograms of target column/s (containing lists) of dataframe dictionary of admixture_times:admixture_types 
    entries 
    '''
    # {{{
    tract_lengths_all_types = []
    labels = []
    for time, typ in dictionary.items():
        tract_lengths = dataframe[(dataframe.admixtype == typ) & (dataframe.t == time)][
            target]  # Subset entries to admixture time and admixture type
        tract_lengths_list = list(tract_lengths)
        tract_lengths_all_sims = [x for y in tract_lengths_list for x in
                                  y]  # putting together all simulations for the subset

        tract_lengths_all_types.append(tract_lengths_all_sims)
        labels.append('time={};type={}'.format(time, typ))

    # plot list of lists of tractlenths ( one for each admixture type) as a side by side histogram
    plt.hist(tract_lengths_all_types, bins=20, density=True, label=labels)
    plt.xlabel('tract lengths')
    plt.ylabel('Frequency')
    plt.title('tract length distributions')
    # }}}
    return ()


def vis_mles(dataframe, target, dictionary):
    '''
    Plot side-by-side histograms of mle's of target column ( should contain list of tract lengths) given dictionary of time:type pairs 
    '''
    # {{{
    mle_list_lists = []
    labels = []

    for t, typ in dictionary.items():
        mle_list_lists.append(get_mle_list_alltracts(dataframe, target, t, typ))
        labels.append('time={};type={}'.format(t, typ))

    plt.hist(mle_list_lists, density=True, label=labels)

    plt.xlabel('MLE sample distribution')
    plt.ylabel('Frequency')
    plt.title('MLE')

    # }}}
    return ()


# Plotting functions (grid)

def plot_row(dataframe, targets, t, typ, compareList=None, mode='mle'):
    '''
    Plots one row of targets side-by-side with an optional compareList given time and admxiture type. Also requires a list of m's to correctmle for t ( finite chromosome correction )
    '''

#{{{ 
    if compareList is None:
        compareList = []

    rowsize = len(targets)
    f, a = plt.subplots(1, rowsize, figsize=(10, 3))

    for idx, ax in enumerate(a):

        column_to_plot = dataframe[targets[idx]][(dataframe.admixtype == typ) & (dataframe.t == t)]
        tract_lengths_list = list(column_to_plot)

        if mode == 'list':
            plot_list = [x for y in tract_lengths_list for x in y]  # putting together all simulations
        else:
            plot_list = tract_lengths_list

        if all(np.isnan(plot_list)):  # Ignore column if plot_list is empty
            ax.plot()
        else:
            ax.hist([plot_list, compareList], bins='auto', alpha=0.7, rwidth=0.85)
            ax.set_title(targets[idx])

    plt.suptitle('typ={}'.format(typ), x=-0.1, y=0.5)
    plt.tight_layout()
#}}}

    return ()


def grid_plot_mles(dataframe, target, t, types, dims):
    """
    plots mle sampling distributions in a grid format given t and a list of all types to plot for
    """
#{{{ 
    mle_list_lists = []
    labels = []

    for typ in types:
        mle_list_lists.append(get_mle_list_alltracts(dataframe, target, t, typ))
        labels.append('time={};type={}'.format(t, typ))

    f, a = plt.subplots(dims[0], dims[1], figsize=(15, 15))
    a = a.ravel()

    for idx, ax in enumerate(a):
        ax.hist(mle_list_lists[idx], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
        ax.set_title(labels[idx])

    plt.tight_layout()
#}}}

    return ()


def grid_plot_tracts(dataframe, target, t, types, dims):
    '''
    plots histogram of tract lengths in a grid format given t and a list of all types to plot for
    '''
#{{{ 
    f, a = plt.subplots(dims[0], dims[1], figsize=(15, 15))
    a = a.ravel()

    for idx, ax in enumerate(a):
        type_list = dataframe[target][dataframe.admixtype == types[idx]]
        tract_lengths_list = list(type_list)
        plot_list = [x for y in tract_lengths_list for x in y]  # putting together all simulations

        ax.hist(plot_list, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
        ax.set_title(types[idx])

    plt.tight_layout()
#}}}

    return ()
