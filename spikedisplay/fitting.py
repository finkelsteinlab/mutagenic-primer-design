import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x, m, b):
    return m*x + b

def fill_in_bins(x, y, unique_bins):
    """
    For each bin number in <unique_bins> that is not found in <x>,
    append that bin number to <x> and 0 to <y>.

    Return modified, x, y
    """
    for unique_bin in unique_bins:
        if unique_bin not in x:
            x = np.append(x, unique_bin)
            y = np.append(y, 0)
            print(f'Bin {unique_bin} not represented')

    return x, y

def fit_bin_count_per_variant(bamdf, variant_index_col='variant_name', **kwargs):
    """
    For every unique value of <variant_index_col> in the <bamdf> dataframe,
    fit a gaussian to read counts vs. bin_number, annotate the result of the
    fit in <countdf>, then aggregate that df to have one row per variant
    and return the aggregate df

    By default, count number of reads per bin per variant and use this "countdf"
    for fitting. You can also pass a premade countdf as a kwarg
    """
    xcol = kwargs.get('xcol', 'bin_number')    
    countdf = kwargs.get('countdf', None)
    normalize_to_total_reads = kwargs.get('normalize_to_total_reads', True)
    # Do this to fill in gaps with zeros for bins
    # not populated at all by each variant
    unique_bins = kwargs.get('unique_bins', [1, 2, 3, 4])
    # Since <countdf> is a groupby df, the index column
    # is the number of reads (rows) in the original 
    # dataframe 
    ycol = kwargs.get('ycol', 'index')
    constraints = kwargs.get('constraints', (0, [1., 6., 6.]))
    filename = kwargs.get('filename', str(bamdf.source_file.iloc[0]).replace('.csv.gz', '_gaussian_fits_df.csv'))
    reads_count_col = kwargs.get('reads_count_col', 'seq_id')
    # Countdf ives number of times each unique variant appears in
    # each bin number
    count_gb_index = [
        'is_wt',
        'mutant_aa',
        'wt_aa_at_mismatch',
        'mismatch_index_in_full_length_wt_seq',
        'source_file',
        'sample_name',
        'library_name',
        'variant_name',
        'sample_library',
        'bin_number'
    ]
    if countdf is None:
        countdf = bamdf.groupby(count_gb_index).count().reset_index()
    variant_index = [variant_index_col]
    countdf.loc[:, 'gaussian_a'] = np.nan
    countdf.loc[:, 'gaussian_center'] = np.nan
    countdf.loc[:, 'gaussian_stdev'] = np.nan
    countdf.loc[:,'total_variant_reads'] = np.nan

    countdf.set_index(variant_index, inplace=True)
    variant_coords = list(countdf.index.unique())
    exceptions = []

    for i, variant_coord in enumerate(variant_coords):

        x = np.array(countdf.loc[variant_coord, xcol])
        y = np.array(countdf.loc[variant_coord, ycol])
        # Append zeroes to y for any missing bins in x
        x, y = fill_in_bins(x, y, unique_bins)
        total_reads = np.sum(countdf.loc[variant_coord, reads_count_col])
        print(f'Fitting variant {i+1} of {len(variant_coords)} with {total_reads} reads')
        if normalize_to_total_reads:
            y_norm = y/np.sum(y)
        else:
            y_norm = y
        try:
            popt, pcov = curve_fit(gaussian, x, y_norm, bounds=constraints)
            a = popt[0]
            center = popt[1]
            sigma = popt[2]
            countdf.loc[variant_coord, 'gaussian_a'] = a
            countdf.loc[variant_coord, 'gaussian_center'] = center
            countdf.loc[variant_coord, 'gaussian_stdev'] = sigma
            countdf.loc[variant_coord, 'total_variant_reads'] = total_reads
            print(f'Fit variant with center {center}')
        except Exception as e:
            print(f'Fit failed with exception {e}')
            exceptions.append(e)
    # Reduce the data down to one row per variant. Note that aggregating by median
    # means that bin_number column will now give a median bin number for that variant
    # Other data, like total_variant_reads, are just the same number multiple times
    # so the value in the pivot table will be that value
    pivot_index = [
        'is_wt',
        'mutant_aa',
        'wt_aa_at_mismatch',
        'mismatch_index_in_full_length_wt_seq',
        'source_file',
        'sample_name',
        'library_name',
        'variant_name',
        'sample_library',
    ]
    fits_df = countdf.pivot_table(index=pivot_index, aggfunc=np.median).reset_index()
    savepath = os.path.join(os.getcwd(), filename)
    fits_df.to_csv(savepath, index=False)
    print(f'Wrote fits df to {savepath}')
    
    return fits_df