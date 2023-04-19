"""reference genome depends on PHYLOP_FILE"""

from contextlib import suppress

import numpy as np
import pyBigWig

from boostdm.globals import PHYLOP_FILE

bw = pyBigWig.open(PHYLOP_FILE)


def get_value(chrom, pos):
    chrom = str(chrom)
    if len(chrom) < 3:
        chrom = 'chr' + str(chrom)
    with suppress(Exception):
        return np.float16(bw.values(chrom, pos - 1, pos)[0])
    return None


def add_feature(df):
    df['PhyloP'] = df.apply(lambda x: get_value(x['chr'], int(x['pos'])), axis=1)
    return df


def func_test():
    """
    output should be:
    -----------------
    Chromosome=19
    Position=1207203
    Phylop=9.53125

    """

    chrom = 'chr10'
    pos = 8055860
    v = get_value(chrom, pos)
    print('Chromosome={0}\nPosition={1}\nPhylop={2}'.format(chrom, str(pos), v))


if __name__ == '__main__':
    func_test()
