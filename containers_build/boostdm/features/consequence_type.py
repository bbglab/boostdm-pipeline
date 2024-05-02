from functools import partial

from boostdm.globals import TABIX_FILE, CONSEQUENCES_LIST, CONSEQUENCES_DICT, AGGREGATION_DICT
from boostdm.vepreader import Tabix



def get_csqn_type(chr_, pos, alt, gene, reader):

    for data in reader.get(chr_, pos, pos):
        
        alt_vep = (data['ALT'] == alt)           # same alternate allele
        mane_vep = (data["MANE_SELECT"] != '-')  # impose mane transcript
        correct_gene = (data["SYMBOL"] == gene)    # skip cases with antisense overlapping genes
        if alt_vep and mane_vep and correct_gene:
            csqn = CONSEQUENCES_LIST[min([CONSEQUENCES_DICT[c] for c in data["CNSQ"].split(',')])]
            return AGGREGATION_DICT.get(csqn, None)
    
    return None


def add_feature(df):

    with Tabix(TABIX_FILE) as reader:
        get_from_reader = partial(get_csqn_type, reader=reader)
        df['csqn_type'] = df.apply(lambda row: get_from_reader(str(row['chr']), int(row['pos']), row['alt'], row['gene']), axis=1)
    return df
