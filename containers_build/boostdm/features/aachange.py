
from functools import partial

import pandas as pd

from boostdm.globals import TABIX_FILE
from boostdm.vepreader import Tabix


def get_aachange(chr_, pos, alt, gene, reader):

    for data in reader.get(chr_, pos, pos):
        alt_vep = (data['ALT'] == alt)
        mane_vep = (data['MANE_SELECT'] != '-') # impose MANE transcript
        correct_gene = (data['SYMBOL'] == gene) # skip cases with antisense overlapping gene (gene is gene_symbol)
        if alt_vep and mane_vep and correct_gene:
            aas = data['AA']  # [11] -> amino-acids involved in change ("I/T")
            aa_pos = data['PROT_POS']  # [10] -> amino-acid position
            if '/' in aas:
                aa_ref, aa_alt = tuple(aas.split('/'))
                return aa_ref + aa_pos + aa_alt
            elif aas == '-':
                return None
            else:
                return aas + aa_pos + aas
    return None


def add_feature(df):

    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_aachange, reader=reader)
        df['aachange'] = df.apply(lambda row: get_from_reader(str(row['chr']), int(row['pos']), row['alt'], row['gene']), axis=1)
    return df


def test():
    """Test function"""
    df = pd.DataFrame({
        'chr': ['14'],
        'pos': [104773083],
        'alt': ['A'],
       'gene': ['AKT1']
    })
    df = add_feature(df)
    print(df)


if __name__ == '__main__':
    test()
