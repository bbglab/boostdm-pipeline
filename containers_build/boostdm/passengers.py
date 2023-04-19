# Usage
# -----


# Imports
# -------

import os
import itertools
from collections import namedtuple
from contextlib import suppress

import bgreference
import numpy as np
import pandas as pd

from boostdm.globals import CANONICAL_TRANSCRIPTS_FILE, GENOME_BUILD


# Utils
# -----

CB = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
TRIPLETS = [p[0] + c + p[1] for c in ['C', 'T'] for p in itertools.product(CB.keys(), repeat=2)]
CHANGES = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}

cds_data = pd.read_csv(CANONICAL_TRANSCRIPTS_FILE,
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

    """Gives index of triplet according to TRIPLET sorting"""
    if triplet[1] not in ['C', 'T']:
        triplet = CB[triplet[2]] + CB[triplet[1]] + CB[triplet[0]]
    return TRIPLETS.index(triplet)


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


def sequence_triplet_index(segments):

    """Takes segments and returns single sequence of triplet indices"""
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


def give_alt_allele(nucleotide, n):

    """Gives the n-th possible change in lexicographic order from nucleotide"""
    return CHANGES[nucleotide][n]


# Indices of reference triplets in standard 96-channel order
# ----------------------------------------------------------

TRIPLET_INDEX_SPECTRA = [triplet_index(triplet) for triplet, alt in mut_key_gen()]


# Probability
# -----------

def probability_vector(sequence, total_mutrate):

    unfold_rates = []
    for ind in sequence:
        if ind == -1:
            unfold_rates += [0, 0, 0]
        else:
            pos = TRIPLET_INDEX_SPECTRA.index(ind)
            for i in range(3):
                unfold_rates.append(total_mutrate[pos + i * 16])
    return np.array(unfold_rates) / sum(unfold_rates)  # TODO: beware of zero-values!


def to_mutation(mut_index, cds, sequence, chrom):

    """
    Return:
        position and index of spectra context
    """

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


# Main: randomize per gene
# ------------------------

def randomize(mutrate, chrom, cds, segments, n_randomizations):

    """
    Return:
        chromosome, Mutation('pos', 'ref_triplet', 'alt')
    """

    sequence = sequence_triplet_index(segments)  # sequence of triplet, given as their indices in TRIPLETS list
    prob = probability_vector(sequence, mutrate)
    try:
        random_draw = np.random.choice(np.arange(len(prob)), size=n_randomizations, p=prob)
    except ValueError:
        prob = 1 / (3 * len(sequence)) * np.ones(3 * len(sequence))  # prob: extended sequence with x3 as many positions
        random_draw = np.random.choice(np.arange(len(prob)), size=n_randomizations, p=prob)

    for mut_index in random_draw:
        with suppress(AssertionError):
            mut = to_mutation(mut_index, cds, sequence, chrom)
            yield chrom, mut


if __name__ == '__main__':

    pass
