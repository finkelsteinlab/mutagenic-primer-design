import numpy as np
import pandas as pd

def flag_normalize(df1, df2, reads_cutoff=20, **kwargs):
    xvar = kwargs.get('xvar', 'gaussian_center_x')
    yvar = kwargs.get('yvar', 'gaussian_center_y')

    df1 = df1[df1.total_variant_reads >= reads_cutoff]
    df2 = df2[df2.total_variant_reads >= reads_cutoff]
    mergedf = df1.merge(df2, on='variant_name', how='left')
    signal_center = mergedf[xvar]
    control_center = mergedf[yvar]

    normalized = signal_center/control_center

    mergedf.loc[:, 'ace2_center_over_flag_center'] = normalized
    
    return mergedf

def wt_normalize(fitsdf, **kwargs):
    
    wt_variant_name = kwargs.get('wt_variant_name', 'None0None')
    colname= kwargs.get('colname', 'gaussian_center')
    colname_ratio = f'{colname}_over_wt'
    colname_log_ratio = f'{colname}_log2_over_wt'

    wt_center_bin = fitsdf.loc[fitsdf.variant_name==wt_variant_name,
                               colname].iloc[0]

    fitsdf.loc[:, colname_ratio] = np.nan
    fitsdf.loc[:, colname_log_ratio] = np.nan

    fitsdf.loc[:, colname_ratio] = fitsdf[colname]/wt_center_bin
    fitsdf.loc[:, colname_log_ratio] = np.log2(fitsdf[colname_ratio])
    
    return fitsdf