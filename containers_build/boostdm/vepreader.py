import tabix


GENOME_SEQUENCE_MAPS = {'chr{}'.format(c): '{}'.format(c) for c in range(1, 23)}
GENOME_SEQUENCE_MAPS.update({'chrX': 'X', '23': 'X', 'chr23': 'X', 'chrY': 'Y', '24': 'Y', 'chr24': 'Y'})
GENOME_SEQUENCE_MAPS.update({'chrM': 'M', 'MT': 'M', 'chrMT': 'M'})


class Tabix:

    def __init__(self, file):
        self.file = file
        self.tb = None
        self.map = GENOME_SEQUENCE_MAPS

    def __enter__(self):
        self.tb = tabix.open(self.file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return True

    def get(self, chromosome, start, stop):
        chr_ = self.map.get(chromosome, chromosome)
        for row in self.tb.query("{}".format(chr_), start, stop):
            yield row

