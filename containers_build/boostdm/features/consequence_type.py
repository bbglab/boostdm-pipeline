from functools import partial

from boostdm.globals import TABIX_FILE, CONSEQUENCES_LIST, CONSEQUENCES_DICT, AGGREGATION_DICT
from boostdm.vepreader import Tabix



def get_csqn_type(chr_, pos, alt, gene, reader):

    for data in reader.get(chr_, pos, pos):
        alt_vep = (data[3] == alt)
        canonical_vep = (data[-4] == 'YES')  # impose canonical transcript
        correct_gene = (data[-7] == gene)  # skip few cases with antisense overlapping gene
        if alt_vep and canonical_vep and correct_gene:
            csqn = CONSEQUENCES_LIST[min([CONSEQUENCES_DICT[i] for i in data[7].split(',')])]
            return AGGREGATION_DICT.get(csqn, None)
    return None


def add_feature(df):

    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_csqn_type, reader=reader)
        df['csqn_type'] = df.apply(lambda row: get_from_reader(str(row['chr']), int(row['pos']), row['alt'], row['gene']), axis=1)
    return df
