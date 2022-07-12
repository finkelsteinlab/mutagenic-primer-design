import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x, m, b):
    return m*x + b

def line_no_intercept(x, m):
    return m*x

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

def undersample_bamdf(bamdf, reads_per_bin=None, use_weights=False, **kwargs):
    """
    Randomly select an equal number of reads from each bin (defaults to lowest number of reads
    per individual bin) and return the undersampled bamdf
    """
    # Aggregate the bamdf to find number of reads per bin
    readsperbindf = bamdf.groupby(by='bin_number').count()
    readsperbindf.loc[:, 'fraction_all_reads_in_bin'] = readsperbindf.seq_id/readsperbindf.seq_id.sum()
    readsperbindf.loc[:, 'bin_sampling_weight'] = 1/readsperbindf.fraction_all_reads_in_bin
    readsperbindf.loc[:, 'frac_bin_sampling_weight'] = readsperbindf.bin_sampling_weight/np.sum(readsperbindf.bin_sampling_weight)
    
    bamdf.loc[:, 'frac_bin_sampling_weight'] = np.nan
    bamdf.loc[:, 'bin_sampling_weight'] = np.nan
    
    for bin_number in list(readsperbindf.index):
        weight = readsperbindf.loc[bin_number, 'bin_sampling_weight']
        frac_weight = readsperbindf.loc[bin_number, 'frac_bin_sampling_weight']
        
        bamdf.loc[bamdf.bin_number==bin_number, 'bin_sampling_weight'] = weight
        bamdf.loc[bamdf.bin_number==bin_number, 'frac_bin_sampling_weight']  = frac_weight
        
    bamdf.loc[:, 'frac_bin_sampling_weight'] = bamdf['bin_sampling_weight']/np.sum(bamdf['bin_sampling_weight'])

    minreadsperbin = readsperbindf.seq_id.min()
    n_bins = len(readsperbindf.index.unique())
    
    if reads_per_bin==None:
        reads_per_bin_to_sample = minreadsperbin

    total_reads_to_sample = n_bins*reads_per_bin_to_sample
    print(f'Sampling {reads_per_bin_to_sample} reads per bin ({total_reads_to_sample} total)')
    sampled_bamdf = bamdf.sample(n=total_reads_to_sample, weights='frac_bin_sampling_weight')
    
    # Sample without using weights
    if not use_weights:
        sampled_dfs = []
        print(f'Sampling bins independently (not whole dataset using weights)')
        for bin_number in list(readsperbindf.index):
            sampled_df = bamdf[bamdf.bin_number==bin_number].sample(reads_per_bin_to_sample)
            sampled_dfs.append(sampled_df)
        sampled_bamdf = pd.concat(sampled_dfs, sort=False)

    return sampled_bamdf

def fit_bin_count_per_variant(bamdf, variant_index_col='variant_name', savedf=True, **kwargs):
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
    if savedf:
        fits_df.to_csv(savepath, index=False)
        print(f'Wrote fits df to {savepath}')
    else:
        print(f'Not saving fits df')
    
    return fits_df

def get_r_squared(y, y_pred):
    """
    Return R-squared as float defined as:

    r_sq = sum of squared model variation from mean / sum of squared data varation from mean

    Note that r-squared is not a valid qualify of fit indicator for non-linear
    regression models. Instead use std. error of the estimate (regression) implemented here
    as get_estimate_std_error
    """
    if (len(y) == len(y_pred)):
        try:
            # Should always work unless y and y_pred aren't np.arrays or
            # pd.DataFrame columns (ie pd.Series)
            y_bar = y.mean()
            r_sq = 1 - (np.power(y - y_pred, 2).sum()) / np.power(y - y_bar, 2).sum()
        except:
            print('Failed to calculate r_sq, y and y_pred were not arrays')
            r_sq = False
    else:
        print("Lengths of y and y_pred not the same")
        r_sq = False

    return r_sq



def normalize_bins(bamdf, cell_counts):
    """
    Return a count df of <bamdf> grouped by variant_name 
    and bin_number with various annotation
    """
    
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
    countdf = bamdf.groupby(count_gb_index).count().reset_index()
    countdf.loc[:, 'bin_enrichment'] = np.nan
    countdf.loc[:, 'fraction_all_reads_in_bin'] = np.nan
    countdf.loc[:, 'fraction_all_cells_in_bin'] = np.nan
    countdf.loc[:, 'total_reads_in_bin'] = np.nan
    countdf.loc[:, 'total_cells_in_bin'] = np.nan

    # Calculate number of reads per bin etc. in the whole data set
    readsperbin = bamdf.groupby(['bin_number']).count()
    readsperbin.loc[:, 'number_of_cells'] = cell_counts
    readsperbin.loc[:, 'nreads_to_ncells'] = readsperbin.seq_id/readsperbin.number_of_cells
    readsperbin.loc[:, 'number_of_cells'] = cell_counts

    total_reads_in_sample = readsperbin['seq_id'].sum()
    total_cells_in_sample = readsperbin['number_of_cells'].sum()
    for bin_number in [1, 2, 3, 4]:

        total_reads_in_bin = readsperbin.loc[bin_number, 'seq_id']
        total_cells_in_bin = readsperbin.loc[bin_number, 'number_of_cells']
        bin_fraction_reads = total_reads_in_bin/total_reads_in_sample
        bin_fraction_cells = total_cells_in_bin/total_cells_in_sample
        # Annotate fraction of reads in bin for this bin
        boolean = countdf.bin_number==bin_number
        countdf.loc[boolean, 'fraction_all_reads_in_bin'] = bin_fraction_reads
        countdf.loc[boolean, 'fraction_all_cells_in_bin'] = bin_fraction_cells
        countdf.loc[boolean, 'total_reads_in_bin'] = total_reads_in_bin
        countdf.loc[boolean, 'total_cells_in_bin'] = total_cells_in_bin

    countdf.loc[:, 'reads_per_cell'] = countdf.total_reads_in_bin/countdf.total_cells_in_bin
    
    return countdf

def get_n_bins_occupied(subdf, warn=False):
    empty_bins = []
    unique_bins = [1, 2, 3, 4]
    for bin_number in unique_bins:
        if bin_number not in subdf.bin_number.values:
            if warn:
                print(f'bin {bin_number} not represented')
            empty_bins.append(bin_number)

    n_bins = len(unique_bins) - len(empty_bins)
    
    return n_bins

def analyze_variant_bin_dists(countdf):

    xcol = 'bin_number'
    ycol = 'seq_id'
    variant_index = 'variant_name'
    # filter out variants that have less than x reads
    countdf = countdf[countdf.seq_id>=0]
    countdf.sort_values(by=xcol, ascending=True, inplace=True)
    if variant_index in countdf.columns:
        countdf.set_index(variant_index, inplace=True)
    variant_coords = list(countdf.index.unique())

    countdf.loc[:, 'fraction_variant_reads_in_bin'] = np.nan
    countdf.loc[:, 'variant_appears_in_n_bins'] = np.nan
    countdf.loc[:, 'max_enriched_bin'] = np.nan
    countdf.loc[:, 'mean_enriched_bin'] = np.nan
    countdf.loc[:, 'bin_enrichment'] = np.nan
    countdf.loc[:, 'bin_enrichment_norm'] = np.nan
    countdf.loc[:, 'gaussian_a'] = np.nan
    countdf.loc[:, 'gaussian_center'] = np.nan
    countdf.loc[:, 'gaussian_stdev'] = np.nan
    countdf.loc[:, 'total_variant_reads'] = np.nan
    countdf.loc[:, 'cell_count_norm_variant_reads_in_bin'] = np.nan
    countdf.loc[:, 'cell_count_norm_variant_fraction_in_bin'] = np.nan
    countdf.loc[:, 'mean_cell_count_norm_bin'] = np.nan
    countdf.loc[:, 'variant_reads_in_bin_over_fraction_all_reads_in_bin'] = np.nan

    variant_count_dfs = []

    for i, variant_coord in enumerate(variant_coords):
        print(f'Analyzing bin distribution of variant {i+1} of {len(variant_coords)}', end='\r')
        # Only quantify if variant has reads in more than one bin
        if not type(countdf.loc[variant_coord, 'bin_number']) == pd.core.series.Series:
            pass
        # If less than two bins, type would be integer
        else:
            # Separate variant slice of countdf into its own dataframe. Makes calculations
            # etc. a little simpler
            subdf = countdf.loc[variant_coord, :].reset_index()
            n_bins = get_n_bins_occupied(subdf)
            subdf.loc[:, 'variant_appears_in_n_bins'] = n_bins
            subdf.loc[:, 'fraction_variant_reads_in_bin'] = subdf.loc[:, ycol]/subdf.loc[:, ycol].sum()
            subdf.loc[:, 'bin_enrichment'] = subdf['fraction_variant_reads_in_bin']/subdf['fraction_all_reads_in_bin']
            subdf.loc[:, 'bin_enrichment_norm'] = subdf['bin_enrichment']/subdf['bin_enrichment'].sum()
            # Weighted mean of bin_enrichment_norm
            mean_enriched_bin = np.sum(subdf['bin_enrichment_norm']*subdf['bin_number'])
            # Discrete bin that is most highly over represented for this variant
            max_enriched_bin = subdf.loc[subdf['bin_enrichment'].idxmax(), 'bin_number']
            max_reads_bin = subdf.loc[subdf[ycol].idxmax(), 'bin_number']
            subdf.loc[:, 'total_variant_reads'] = subdf[ycol].sum()
            subdf.loc[:, 'mean_enriched_bin']  = mean_enriched_bin
            subdf.loc[:, 'max_enriched_bin'] = max_enriched_bin
            subdf.loc[:, 'max_reads_bin'] = max_reads_bin
            subdf.loc[:, 'mean_bin'] = np.sum(subdf['fraction_variant_reads_in_bin']*subdf['bin_number'])
            # Normalize reads per bin to fraction all reads in bin. This will amplify low read bins
            subdf.loc[:, 'variant_reads_in_bin_over_fraction_all_reads_in_bin'] = subdf.loc[:, ycol]/subdf.loc[:, 'fraction_all_reads_in_bin']
            subdf.loc[:, 'variant_reads_in_bin_amplified_norm'] = subdf['variant_reads_in_bin_over_fraction_all_reads_in_bin']/np.sum(subdf['variant_reads_in_bin_over_fraction_all_reads_in_bin'])
            subdf.loc[:, 'weighted_mean_amplified_bin'] = np.sum(subdf['variant_reads_in_bin_amplified_norm']*subdf['bin_number'])
            max_amplified_bin = subdf.loc[subdf['weighted_mean_amplified_bin'].idxmax(), 'bin_number']
            subdf.loc[:, 'max_amplified_bin'] = max_amplified_bin
            # Another potential normalization scheme based on Matouschek lab procedure
            subdf.loc[:, 'cell_count_norm_variant_reads_in_bin'] = (subdf['fraction_all_cells_in_bin']/subdf['total_reads_in_bin'])*subdf[ycol]
            subdf.loc[:, 'cell_count_norm_variant_fraction_in_bin'] = subdf['cell_count_norm_variant_reads_in_bin']/np.sum(subdf['cell_count_norm_variant_reads_in_bin'])
            subdf.loc[:, 'mean_cell_count_norm_bin'] = np.sum(subdf['cell_count_norm_variant_fraction_in_bin']*subdf[xcol])
            # Gaussian fitting that is of dubious use
            constraints = (0, [1., 4., 4.])
            try:
                popt, pcov = curve_fit(gaussian, subdf.bin_number.values, subdf.bin_enrichment_norm.values, bounds=constraints)
                a = popt[0]
                center = popt[1]
                sigma = popt[2]
                
                subdf.loc[:, 'gaussian_a'] = a
                subdf.loc[:, 'gaussian_center'] = center
                subdf.loc[:, 'gaussian_sigma'] = sigma
            except:
                pass
            variant_count_dfs.append(subdf)
    countdf = pd.concat(variant_count_dfs, ignore_index=True)
    return variant_count_dfs, countdf