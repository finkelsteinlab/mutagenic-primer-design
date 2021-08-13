import os
import datetime
from itertools import compress

import pandas as pd
import numpy as np
# BioPython imports
import Bio
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp as mt

import Levenshtein
# spikedisplay imports
from spikedisplay import constants

def clean_NNN_primers(output_path):
    """
    Return the .txt CSV found at *output* path with
    blank spaces removed and reverse primers dropped
    """
    primer_df = pd.read_csv(output_path)
    primer_df.reset_index(inplace=True)
    primer_df.columns = ['primer_name', 'sequence_nnn']
    # Drop nan values (empty spaces between plates)
    nan_mask = ~primer_df.sequence_nnn.isna()
    primer_df = primer_df.loc[nan_mask, :]
    # Drop reverse primers
    rev_mask = ['rev' not in name for name in primer_df.primer_name]
    primer_df = primer_df.loc[rev_mask, :]
    primer_df.index = np.arange(0, len(primer_df), 1)

    return primer_df

def annotate_NNN_primers(primer_df, sequence_file_path):
    """
    Take primer_df, add a column for NN melt temp, WT codon,
    and WT aminio acid at each codon
    """
    # Read in the sequence file and select the open reading frame
    # (ORF) recognized as all upper case characters in the file
    with open(sequence_file_path) as file:
        seq = file.read()
        seq_list = list(seq)
        bools = [string.isupper() for string in seq_list]
        ORF_list = list(compress(seq_list, bools))
        ORF = ''.join(ORF_list)

    # Split the open reading frame (ORF) into a list of codons
    codons = [ORF[i:i+3] for i in range(0, len(ORF), 3)]
    # Add original codon that NNN replaced to primer dataframe (*primer_df*)
    primer_df.loc[:, 'codon_to_replace'] = codons
    # Add amino acid abbreviation to primer dataframe
    amino_acids = [str(Seq(codon).translate()) for codon in codons]
    primer_df.loc[:, 'amino_acid'] = amino_acidss
    # Add melt temp to primer dataframe
    primers = [Seq(p) for p in list(primer_df.sequence_nnn)]
    mts = [mt.Tm_NN(p) for p in primers]
    primer_df.loc[:, 'Tm_NN'] = mts

    return primer_df

def evaluate_aa_codons(primer_df, primer_idx):
    """
    *primer_df*: a dataframe with columns 'primer_name',
    'sequence_nnn', 'codon_to_replace', 'amino_acid'] generated
    with 
    
    *primer_idx*: dataframe index of the primer for which we want to 
    replace "NNN" with an optimal codon for each of the 20 amino acids
    
    Return a list of dataframes, one for amino acid of 20. Each contains
    actual primer melt temp, Levenshtein distance from original codon, codon
    usage in humans in terms of fraction, etc
    """
    print(f'Finding alternative codons for primer number {primer_idx}')
    primer_seq_NNN = primer_df.loc[primer_idx, 'sequence_nnn']
    primer_Tm_NN = primer_df.loc[primer_idx, 'Tm_NN']
    codon_table = pd.read_csv(constants.human_codon_table_path)
    # Don't need stop codons in this library so will drop all
    # stop codons from the *codon_table*
    codon_table = codon_table[codon_table.amino_acid!='*']
    # Iterating through primers in the generated NNN library
    codon_to_replace = primer_df.codon_to_replace[primer_idx]
    orig_aa = primer_df.amino_acid[primer_idx]
    # Iterate through all amino acids in the *codon_table*, find the
    # Codon for each with the highest distance from *codon_to_replace*
    # and the highest codon usage
    codons_dfs = []
    for aa in list(codon_table.amino_acid.unique()):
        print(f'Evaluating codons for amino acid {aa}')
        cols = ['triplet', 'fraction']

        if len(codon_table[codon_table.amino_acid==aa]) == 0:
            print(f'No codons found for amino acid {aa}')
            codons_df = None
        elif len(codon_table[codon_table.amino_acid==aa]) == 1:
            print(f'Only one codon for amino acid {aa}')
            codons_df = pd.DataFrame(index=[0], columns=cols)
            codons_df.loc[:, 'amino_acid'] = aa
            new_codon = codon_table.set_index(['amino_acid']).loc[aa, 'triplet']
            codons_df.loc[:, 'fraction'] = codon_table.set_index(['amino_acid']).loc[aa, 'fraction']
            codons_df.loc[:, 'triplet'] = new_codon
            dists = Levenshtein.distance(codon_to_replace, new_codon)
        else:        
            codons_df = codon_table.set_index(['amino_acid']).loc[aa, cols].reset_index()
            # Add Levenshtein distance from original codon to codons_df for this amino acid
            dists = [Levenshtein.distance(codon_to_replace, codon) for codon in codons_df.triplet]
        # Add Levenshtein distance to dataframe
        codons_df.loc[:, 'levenshtein_dist'] = dists
        # Label the dataframe with the original codon and original amino acid we're trying
        # to mutate from
        codons_df.loc[:, 'original_codon'] = codon_to_replace
        codons_df.loc[:, 'original_amino_acid'] = orig_aa
        codons_df.loc[:, 'primer_number'] = primer_idx
        codons_df.loc[:, 'primer_nnn'] = primer_seq_NNN
        codons_df.loc[:, 'primer_Tm_NN'] = primer_Tm_NN
        codons_df.loc[:, 'primer_orig'] = primer_seq_NNN.replace('NNN', codon_to_replace)
        # Annotate what each primer sequence would be for each candidate codon for this amino acid
        # then calculate actual melt anneal temp to original sequence
        new_primer_seqs = [primer_seq_NNN.replace('NNN', seq) for seq in codons_df.triplet]
        codons_df.loc[:, 'primer_new'] = new_primer_seqs
        codons_df.loc[:, 'Tm_real'] = np.nan
        for idx in list(codons_df.index):
            orig_seq = codons_df.loc[idx, 'primer_orig']
            new_seq = codons_df.loc[idx, 'primer_new']
            # if orig_seq == new_seq:
            #     print(f'WARNING: generated new sequence same as original codon')
            # else:
            #     print(f'Successfully generated new codon')
            Tm = mt.Tm_GC(orig_seq, new_seq)
            codons_df.loc[:, 'Tm_real'] = Tm
        # Add codons_df for this amino acid to the list of dfs for this primer position
        # to be returned to the user
        codons_dfs.append(codons_df)
    n_amino_acids_evaluated = len(codons_dfs)
    print(f'\nEvaluated codons of {n_amino_acids_evaluated} amino acids for primer number {primer_idx}')
    
    return codons_dfs

def score_and_select_codon_dfs(primer_codon_dfs, **kwargs):
    """
    Default scoring method is to just add up Levenshtein distance
    from original sequence and codon usage in terms of fraction.
    This naturally weights levenshtein distance (0-3) more highly 
    than codon usage (0-1)
    """
    # Choose which ranked codon to select. They are
    # ranked in order from highest score to lowest
    rank_to_select = kwargs.get('rank_to_select', 0)
    return_scored_df = kwargs.get('return_scored_df', False)
    usage_weight = kwargs.get('usage_weight', 1) 
    mismatch_weight = kwargs.get('mismatch_weight', 1)
    dfs = primer_codon_dfs
    selected_df = pd.DataFrame()
    selected_primers = []
    scored_dfs = []
    for df in dfs:
        df.loc[:, 'score'] = df.fraction*usage_weight + df.levenshtein_dist*mismatch_weight
        sorted_df = df.sort_values(by='score', ascending=False, ignore_index=True)
        scored_dfs.append(sorted_df)
        if len(sorted_df) == 1:
            # There's only one codon for this amino acid so we
            # only have one choice
            selected_primer = sorted_df.loc[0, 'primer_new']
        else:
            selected_primer = sorted_df.loc[rank_to_select, 'primer_new']
        selected_primers.append(selected_primer)
    print(f'Selected {len(selected_primers)} primers for this site with rank {rank_to_select}')
    if return_scored_df:
        return (selected_primers, scored_dfs)
    else:
        return selected_primers

def generate_library(primer_df, rank_to_select, keep_wt=True,
                     return_dfs=False, writepath=None, **kwargs):
    """
    Take the primer_df, which contains a list of NNN primers generated by
    J. Bloom's NNN codon tiling script, generate lists of primers to mutate
    each codon in the primer_df to every amino acid, select the ones with
    highest (or other *rank_to_select*) additive score of levenshtein distance
    from original codon and codon usage in terms of fraction.

    Return the selected primers df that can be used in final library

    If *return_unselected_df*, also return the raw df of all primers found
    for every amino acid at every position before selection
    """
    usage_weight = kwargs.get('usage_weight', 1)
    mismatch_weight = kwargs.get('mismatch_weight', 1)
    timestamp = datetime.datetime.today().timestamp()
    # Set a name for the primer library .csv to be saved at.
    if writepath == None:
        library_name = primer_df.primer_name.iloc[0].split('-')[0]
        weights = f'usage-weight-{usage_weight}_mismatch-weight-{mismatch_weight}'
        filename = f'{library_name}_rank-{rank_to_select}_{weights}_oligo_library_{timestamp}.csv'
        writepath = os.path.join(constants.source_path, filename)
    all_aa_codons_dfs = []
    selected_primer_dfs = []
    # Set keyword arguments for how to score and select
    # the final library of primers
    score_kwargs = {'rank_to_select': rank_to_select,
                    'return_scored_df': True,
                    'usage_weight': usage_weight,
                    'mismatch_weight': mismatch_weight,}
    for primer_idx in list(primer_df.index):

        amino_acid_dfs = evaluate_aa_codons(primer_df, primer_idx)
        primers, scored_dfs = score_and_select_codon_dfs(amino_acid_dfs, **score_kwargs)
        # Combine all data collected on each primer and mutants to return to user if desired
        all_aa_codons_df = pd.concat(scored_dfs, ignore_index=True)
        all_aa_codons_dfs.append(all_aa_codons_df)
        # Selected top ranked primers output by <score_and_select_codon_dfs() above>
        mask = all_aa_codons_df.primer_new.isin(primers)
        selected_primer_df = all_aa_codons_df[mask].reset_index()
        selected_primer_dfs.append(selected_primer_df)
    # Combine unselected primers at each codon into one dataframe
    unselected_primers_df = pd.concat(all_aa_codons_dfs, ignore_index=True)
    # Combine the selected primers at each codon into one dataframe
    selected_primers_df = pd.concat(selected_primer_dfs, ignore_index=True).reset_index()
    if not keep_wt:
        print('Dropping codons for WT amino acids')
        mask = selected_primers_df.amino_acid != selected_primers_df.original_amino_acid
        selected_primers_df = selected_primers_df[mask]
        selected_primers_df.index = np.arange(0, len(selected_primers_df), 1)
        writepath = writepath.replace('_oligo_library_', '_oligo_library_no-WT_')
    else:
        print('Keeping codons for WT amino acid')
    # Annotate the final library dataframe with information about how
    # it was created (weights for mismatch and codon usage etc.)
    selected_primers_df.loc[:, 'codon_usage_weight'] = usage_weight
    selected_primers_df.loc[:, 'mismatch_weight'] = mismatch_weight
    selected_primers_df.loc[:, 'selected_primer_rank'] = rank_to_select
    selected_primers_df.to_csv(writepath, index=False)
    print(f'Saved primer library at {writepath}')

    if return_dfs:
        return (selected_primers_df, unselected_primers_df)