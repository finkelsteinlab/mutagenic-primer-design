import os

import numpy as np
import pandas as pd
# Plotting tools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

plt_params_dict = {'font.sans-serif': 'Arial',
                   'mathtext.default': 'regular'}

def set_styles(plt, matplotlib):
    """
    Set fonts and default style
    """

    try:
        plt.style.use('default')        
        for param, value in plt_params_dict.items():
            print(param, value)
            plt.rcParams.update({param: value})
    except:
        print("""
            Before running set_styles(), you must:

            import matplotlib.pyplot as plt
            import matplotlib
            """)

def legend_outside():
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

def human_format(
    ax,
    axname='x',
    divide_by=1000,
    suffix='k',
    **kwargs
):
    """
    Make axes with very large numbers human readable

    Returns nothing
    """

    def func(x, pos):
        'The two args are the value and tick position'
        new_number = int(x/divide_by)
        new_text = f'{new_number}{suffix}'
        return new_text

    formatter = FuncFormatter(func)

    if axname=='x':
        ax.xaxis.set_major_formatter(formatter)
    elif axname == 'y':
        ax.yaxis.set_major_formatter(formatter)
    else:
        print(f'Please specify x or y axis')

def format_ticks(ax, **kwargs):
    """
    take <ax>, turn on minor ticks and set all ticks
    to inner direction. Can pass <xminorspace> and
    <yminorspace> as keyword args. Otherwise defaults
    to a minor tick every 1/5 of space between major 
    ticks.

    Return nothing
    """
    format_x_ticks = kwargs.get('format_x_ticks', True)
    format_y_ticks = kwargs.get('format_y_ticks', True)
    tickdirection = kwargs.get('tickdirection', 'in')
    # Default to minor tick every 1/5 of space between
    # major ticks
    try:
        xmajorspace = np.diff(ax.get_xticks())[0]
    except:
        xmajorspace = np.nan
    ymajorspace = np.diff(ax.get_yticks())[0]
    xminorspace = xmajorspace/5
    yminorspace = ymajorspace/5

    yminorspace = kwargs.get('yminorspace', yminorspace)
    xminorspace = kwargs.get('xminorspace', xminorspace)
    
    if format_x_ticks:
        ax.xaxis.set_minor_locator(MultipleLocator(xminorspace))

    if format_y_ticks:
        ax.yaxis.set_minor_locator(MultipleLocator(yminorspace))

    for ticktype in ['minor', 'major']:
        ax.tick_params(axis="y", which=ticktype, direction=tickdirection)
        ax.tick_params(axis="x", which=ticktype, direction=tickdirection)

def remove_spines(ax):
    """
    Hide top and right spines on <ax>
    """
    for spine in [ax.spines[key] for key in ['top', 'right']]:
        spine.set_visible(False)



def plot_library_correlations(plasmiddf, libdf, **kwargs):
    
    filetype = kwargs.get('filetype', 'png')

    fig = plt.figure(figsize=(9, 3))
    fig.set_dpi(200)

    ylim = (0, 1000)
    xlim = (0, 1000)

    scatterkwargs = {
        'edgecolor':'black',
        'facecolor':'white',
        's':3,
        'alpha': 0.2
    }


    i=1
    plasmid_lib_name = plasmiddf.sample_library.unique()[0]
    filename = f'{plasmid_lib_name}_variant_coverage_correlations.{filetype}'

    sample_libraries = [0, 1, 2]
    for lib in libdf.sample_library.unique():
        if 'Integrated' in lib:
            sample_libraries[0] = lib
        elif 'FLAG' in lib:
            sample_libraries[1] = lib
        elif 'ACE2' in lib:
            sample_libraries[2] = lib
    for sample_library in sample_libraries:
        if not sample_library == plasmid_lib_name:
            subdf = libdf[libdf.sample_library==sample_library]
            mergedf = plasmiddf.merge(subdf, how='left', on='variant_name')
            ax = fig.add_subplot(1, 3, i)

            x = mergedf.seq_id_x
            y = mergedf.seq_id_y
            ax.scatter(x, y, **scatterkwargs)
            remove_spines(ax)

            xlabel = f'Reads per variant in {plasmid_lib_name}'
            ylabel = f'Reads per variant in {sample_library}'
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax.set_yscale('log')
            ax.set_xscale('log')

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)    
            i += 1

    plt.tight_layout()
    savepath = os.path.join(os.getcwd(), filename)
    fig.savefig(savepath)
    print(f'Saved figure at {savepath}')
    return ax