import os
import tarfile
import json

import tqdm


class MutrateReader:
 
    def __init__(self, input_file):

        self.input_file = input_file
        self.genes = self._get_genes()

    def _get_genes(self):

        with tarfile.open(self.input_file, 'r|gz') as tf:
            members = tf.getmembers()
            genes_dict = {}
            for m in members:
                info = m.get_info()
                if info['name'].endswith('.json'):
                    gene = os.path.basename(info['name']).split('.')[0]
                    genes_dict[gene] = info['name']
            return genes_dict

    def load(self, symbol):

        with tarfile.open(self.input_file, 'r') as tf:
            try:
                res = json.load(tf.extractfile(self.genes[symbol]))
            except KeyError as err:
                print(err)
                return None
            except OSError as err:
                print('TarFile is closed')
                return None
            return res[symbol]


if __name__ == '__main__':

    """testing MutrateReader"""

    fn = '/workspace/projects/intogen_2017/runs/20200102/mutrate/CBIOP_WXS_ACC_2019.tar.gz'
    drivers = ['AKT1']

    dict_genes = {}
    mr = MutrateReader(fn)
    print("Reading mutrates for driver genes...")
    for gene in tqdm.tqdm(drivers):
        dict_genes[gene] = mr.load(gene)
    print(dict_genes)
