import os
import click

import pandas as pd

# import bgvep
# from bgvep import get
# from bgvep.readers import Tabix

# TODO: not using bgvep because of:
# 1) bgdata.errors.TagNotFound: Tag master for package vep/wgs_tabix/hg38_92 not found
# 2) exon number is not provided as VEP field
# is there any possible workaround?

from boostdm.globals import TABIX_FILE
from boostdm.vepreader import Tabix
# genome build and VEP version are dictated by the VEP TABIX_FILE

from boostdm.globals import CANONICAL_TRANSCRIPTS_FILE

# genome build and VEP version
genome_build = os.environ['GENOME_BUILD']
vep_version = os.environ['VEP_VERSION']

# genomic region comprised
regions_fn = CANONICAL_TRANSCRIPTS_FILE
df_regions = pd.read_csv(regions_fn, sep='\t', low_memory=False)

# drivers summary table
drivers_fn = os.environ['DRIVERS_PATH']
df_drivers = pd.read_csv(drivers_fn, sep='\t', low_memory=False)
drivers = df_drivers['SYMBOL'].unique()

# VEP headers

"""
headers = ['Chromosome', 'Position', 'Reference', 'Alternate', 'Gene', 'Feature',
           'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position',
           'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
           'Impact','Distance', 'Strand', 'Flags', 'Symbol', 'Symbol source',
           'HGNC_ID', 'Canonical', 'ENSP']
"""

tmp_headers = list(range(25))


@click.command()
def run():

    for driver in drivers:
        regions = df_regions[df_regions['SYMBOL'] == driver]
        regions.sort_values(by=['GENE_ID', 'START'], inplace=True)

        mutations = []
        for i, r in regions.iterrows():
            start = r['START']
            end = r['END']
            chr_ = str(r['CHROMOSOME'])
            with Tabix(TABIX_FILE) as reader:
                for data in reader.get('chr' + chr_, int(start), int(end)):
                    if data[21] == 'YES':
                        # only if the segment maps to the canonical transcript
                        mutations.append([x for x in data])

        df = pd.DataFrame(mutations, columns=tmp_headers)
        fn = os.path.join(driver + ".vep.gz")
        df.to_csv(fn, index=False, compression='gzip', sep='\t')


@click.command()
def test():
    for driver_gene in drivers:
        df = pd.DataFrame()
        fn = os.path.join(driver_gene + ".vep.gz")
        df.to_csv(fn, index=False, compression='gzip', sep='\t')


if __name__ == '__main__':
    run()
