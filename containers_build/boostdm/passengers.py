# Usage
# -----


# Imports
# -------

import os
import itertools
from itertools import product
from collections import namedtuple
from contextlib import suppress

import bgreference
import numpy as np
import pandas as pd

from boostdm.globals import MANE_TRANSCRIPTS_FILE, GENOME_BUILD


# Utils
# -----

subs = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]  
flanks = [''.join(z) for z in product('ACGT', repeat=2)]  
contexts_tuple_format = sorted([(s, f) for s, f in product(subs, flanks)], key=lambda x: (x[0], x[1]))  
contexts_mutrate_format = [f[0] + s[0] + f[1] + '>' + s[1] for s, f in contexts_tuple_format]



CB = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
TRIPLETS = [p[0] + c + p[1] for c in ['C', 'T'] for p in itertools.product(CB.keys(), repeat=2)]
CHANGES = {'C': ['A', 'G', 'T'], 'T': ['A', 'C', 'G']}

cds_data = pd.read_csv(MANE_TRANSCRIPTS_FILE,
                       sep='\t', header=None, compression='gzip', low_memory=False)
cds_data = cds_data[[0, 1, 2, 3, 6]].copy()
cds_data.columns = ['chr', 'start', 'end', 'strand', 'gene']


# Retrieve exons
# --------------

def retrieve_exons(gene):

    """
    Returns
        chromosome
        CDS of the gene: list of genomic positions of the CDS
        list of exon sequences (): each exon sequence has 1-bp offset at the flanks
    """

    df = cds_data[cds_data['gene'] == gene].copy()
    exons = []
    cds = []

    if GENOME_BUILD == 'hg38':
        func = bgreference.hg38
    elif GENOME_BUILD == 'hg19':
        func = bgreference.hg19

    for ind, row in df.iterrows():
        exons.append(func(row['chr'], int(row['start']) - 1, size=int(row['end']) - int(row['start']) + 3))
        cds += list(range(int(row['start']), int(row['end']) + 1))
    return row['chr'], cds, exons


# Triplet Utils
# -------------

def triplet_index(triplet):

    """Gives index of triplet according to TRIPLETS sorting"""
    if triplet[1] not in ['C', 'T']:
        triplet = CB[triplet[2]] + CB[triplet[1]] + CB[triplet[0]]
    return TRIPLETS.index(triplet)

"""
def reverse_complement(triplet):

    return CB[triplet[2]] + CB[triplet[1]] + CB[triplet[0]]
"""


def reverse_complement(triplet):

    return CB[triplet[2]] + CB[triplet[1]] + CB[triplet[0]]


def mut_key_gen():

    for ref in ['C', 'T']:
        for alt in CB.keys():
            if ref == alt:
                continue
            else:
                for p in itertools.product(CB.keys(), repeat=2):
                    yield p[0] + ref + p[1], alt


"""
def sequence_triplet_index(segments):

    # Takes segments and returns single sequence of triplet indices
    sequence = []
    for seg in segments:
        for i, c in enumerate(seg[1: -1]):
            triplet = seg[i: i+3]
            tset = set(triplet)
            if not tset.issubset(set('ACGT')):
                sequence.append(-1)  # index = -1 implies that mutrate not defined,
                                     # then we impute in the function probability_vector
            else:
                sequence.append(triplet_index(triplet))
    return sequence
"""


def get_triplet_sequence(segments):

    # Takes segments and returns single sequence of triplet indices
    sequence = []
    for seg in segments:
        for i, c in enumerate(seg[1: -1]):
            triplet = seg[i: i+3]
            tset = set(triplet)
            if not tset.issubset(set('ACGT')):
                sequence.append(None)  # None implies that the mutrate is not defined,
                                       # then we impute in the function probability_vector
            else:
                sequence.append(triplet)
    return sequence


"""
def give_alt_allele(nucleotide, n):

    # Gives the n-th possible change in lexicographic order from nucleotide
    return CHANGES[nucleotide][n]
"""
    

"""
# Indices of reference triplets in standard 96-channel order
# ----------------------------------------------------------

# TRIPLET_INDEX_SPECTRA = [triplet_index(triplet) for triplet, alt in mut_key_gen()]
"""


# Probability
# -----------

def probability_vector(sequence, mutrate):

    """
    mutrate is a dictionary with context keys and numerical values:
    https://bbglab.github.io/bbgwiki/Tools/Signature_tools/TrinucleotideOrdering/
    """

    unfold_rates = []
    for triplet in sequence:
        if triplet is None:
            unfold_rates += [0, 0, 0]
        else:
            if triplet[1] not in {'C', 'T'}:
                triplet = reverse_complement(triplet)
            ref = triplet[1]
            for i, alt in enumerate(CHANGES[ref]):
                unfold_rates.append(mutrate[triplet + '>' + alt])
    return np.array(unfold_rates) / sum(unfold_rates)


"""
def to_mutation(mut_index, cds, sequence, chrom):

    # Return:
    #     position and index of spectra context
    
    index_in_sequence = mut_index // 3
    change = mut_index % 3
    pos = cds[index_in_sequence]

    if GENOME_BUILD == 'hg38':
        func = bgreference.hg38
    elif GENOME_BUILD == 'hg19':
        func = bgreference.hg19

    ref_triplet = func(chrom, pos - 1, size=3)
    seq_triplet = TRIPLETS[sequence[index_in_sequence]]
    assert((ref_triplet == seq_triplet) or (ref_triplet == reverse_complement(seq_triplet)))
    alt = give_alt_allele(ref_triplet[1], change)
    Mutation = namedtuple('Mutation', ['pos', 'ref_triplet', 'alt'])
    return Mutation(pos, ref_triplet, alt)
"""


def to_mutation(mut_index, cds, sequence):

    # Return:
    #     position, reference triplet and alternate allele of the mutation
    
    index_in_sequence = mut_index // 3
    change_index = mut_index % 3
    pos = cds[index_in_sequence]
    ref_triplet = sequence[index_in_sequence]

    if ref_triplet[1] not in {'C', 'T'}:
        triplet_pyr = reverse_complement(ref_triplet)
        alt_pyr = CHANGES[triplet_pyr[1]][change_index]
        alt = CB[alt_pyr]
    else:
        alt = CHANGES[ref_triplet[1]][change_index]
    
    Mutation = namedtuple('Mutation', ['pos', 'ref_triplet', 'alt'])
    return Mutation(pos, ref_triplet, alt)


# Main: randomize per gene
# ------------------------

def randomize(mutrate, chrom, cds, segments, n_randomizations):

    """
    Return:
        chromosome, Mutation('pos', 'ref_triplet', 'alt')
    """

    sequence = get_triplet_sequence(segments)  # sequence of triplets
    prob = probability_vector(sequence, mutrate)
    try:
        random_draw = np.random.choice(np.arange(len(prob)), size=n_randomizations, p=prob, replace=True)
    except ValueError:
        prob = 1 / (len(prob)) * np.ones_like(prob)  # prob: extended sequence with x3 as many positions
                                                     # uniform probability distribution
        random_draw = np.random.choice(np.arange(len(prob)), size=n_randomizations, p=prob, replace=True)

    for mut_index in random_draw:
        with suppress(AssertionError):
            mut = to_mutation(mut_index, cds, sequence)
            yield chrom, mut


if __name__ == '__main__':

    pass
