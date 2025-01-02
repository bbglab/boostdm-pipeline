from functools import partial

import pandas as pd

from boostdm.globals import TABIX_FILE
from boostdm.vepreader import Tabix


def get_nmd(chr_, pos, alt, gene, reader):

    for data in reader.get(chr_, pos, pos):
        alt_vep = (data["ALT"] == alt)
        mane_vep = (data["MANE_SELECT"] != '-') # impose mane transcript
        correct_gene = (data["SYMBOL"] == gene)  # skip cases with antisense overlapping gene
        if alt_vep and mane_vep and correct_gene:
            if data["NMD_SKIPPING"] == '-':
                return 0
            elif data["NMD_SKIPPING"] == 'NMD_escaping_variant':
                return 1
    return 0


def add_feature(df):

    df = df.copy()
    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_nmd, reader=reader)
        df['nmd'] = df.apply(lambda row: get_from_reader(row['chr'], row['pos'], row['alt'], row['gene']), axis=1)
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
