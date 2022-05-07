import os

import pandas as pd

def package_path(**kwargs):
    """
    Return the path to the local installation
    of spikedisplay
    """
    import spikedisplay
    package_path = spikedisplay.__file__
    package_path = os.path.dirname(package_path)
    return package_path

def source_path(**kwargs):
    """
    Return the path to the local installation
    source of spikedisplay. (Includes 'data', 'notebooks',
    'bin', 'docs', etc.)
    """
    import spikedisplay
    dirname = os.path.dirname(spikedisplay.__file__)
    source_path = os.path.dirname(dirname)
    return source_path
    
package_path = package_path()
source_path = source_path()

bin_dir = os.path.join(source_path, 'bin')
human_codon_table_path = os.path.join(package_path, 'Human_Codon_Usage.csv')
rbd_barcode_index = pd.read_csv('/root/projects/dms/data/Barcodes_RBD_NGS.csv')

cut_sites = {'BsmBI': 'CGTCTC',
             'BsaI': 'GGTCTC',
             'AaRI': 'CACCTGC',
             'NotI': 'GCGGCCGC'}

class Patterns(object):

    def __init__(self):

        self.barcode = self.get_fwd_rev_barcodes()

    def get_fwd_rev_barcodes(self):
        """
        Return the pattern used in a regex search
        to identify master index files
        """
        # Has string 'master_index' anyhwere in string,
        # ends with '.csv'
        bases = '|'.join(['G', 'C', 'A', 'T'])
        pattern = ''
        for i in range(24):
            pattern = pattern + f'[{bases}]'
        pair_pattern = f'({pattern})_({pattern})'
        return pair_pattern

patterns = Patterns()