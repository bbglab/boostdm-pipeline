from functools import partial

import pandas as pd

from boostdm.globals import TABIX_FILE
from boostdm.vepreader import Tabix


def nmd_rule(exon, total_exons):

    """
    mutation is in an exon
    if first or last exon then return 1
    otherwise return 0
    """

    if exon == 0:
        nmd = 0
    elif (exon == 1) or (exon == total_exons):
        nmd = 1
    else:
        nmd = 0
    return nmd


def get_exon(chr_, pos, alt,gene, reader):

    for data in reader.get(chr_, pos, pos):
        alt_vep = (data["ALT"] == alt)
        mane_vep = (data["MANE_SELECT"] != '-') # impose mane transcript
        correct_gene = (data["SYMBOL"] == gene)  # skip cases with antisense overlapping gene
        if alt_vep and mane_vep and correct_gene:
            exons = data["EXON"]
            if '/' in exons:
                exon, total_exons = tuple(exons.split('/'))
            else:
                exon, total_exons = 0, 0
            return nmd_rule(exon, total_exons)
    return 0


def add_feature(df):

    df = df.copy()
    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_exon, reader=reader)
        df['nmd'] = df.apply(lambda row: get_from_reader(str(row['chr']),
                                                         int(row['pos']),
                                                         row['alt'],
                                                         row['gene']), axis=1)
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
