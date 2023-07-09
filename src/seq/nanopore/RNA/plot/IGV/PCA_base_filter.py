#!/usr/bin/env python3
# Filters bases, based on PCA results.

import math
import os
import pdb
import re

import numpy as np
import pandas

import Bio.Align
import pysam

PCA_dir = 'PC_scores/'

# read in the clone sequence
# (previously, this was only the clone sequence; now reading it from chr19)
clone_fasta = pysam.FastaFile('../../../../../../chr19_AANCR_clone_only_230411.fa')
clone_seq = clone_fasta.fetch('chr19', 45406984, 45408892)
clone_seq = clone_seq.upper().replace('T', 'U')

def maha_2d_dist_cutoff(p):
    """Euclidean 2-D distance cutoff.

    Based on https://en.wikipedia.org/wiki/Mahalanobis_distance
    """
    return math.sqrt(-2. * math.log(1.-p))

p = 0.9
print(f'2-D cutoff for p = {p} = {maha_2d_dist_cutoff(p)}')

def alignment_blocks_to_positions(aligned_blocks):
    """Converts from a list of blocks to a list of positions."""
    # first, expand each block into a list of positions
    pos1 = [list(range(b[0], b[1])) for b in aligned_blocks]
    # then, concatenate these lists
    return sum(pos1, [])

def get_pos_in_clone_seq(clone_seq, seq):
    """Gets positions in the clone sequence, by aligning.

    clone_seq: the clone sequence
    seq: the bases to be matched up
    Returns: a Numpy array, giving the aligned (1-based) indices
        in clone_seq
    """
    # get the best alignment
    aligner = Bio.Align.PairwiseAligner()
    # aligner.mode = 'local'
    aligner.mismatch_score = 0.5
    # aligner.gap_score = -1
    alignments = aligner.align(clone_seq, seq)
    best_alignment = next(alignments)
    (target_pos, query_pos) = [alignment_blocks_to_positions(b)
        for b in best_alignment.aligned]
    # note that this should return None for any of these which are missing
    query_to_target = dict(zip(target_pos, query_pos))
    # get positions for all sites in the query sequence
    def get_pos(i):
        if i in query_to_target:
            return query_to_target[i]
        else:
            return None
    pos = [get_pos(i) for i in range(len(seq))]
    # ??? FIXME add check that the sequences match up?
    return pos

def parse_PCA(filename):
    """Parses one PCA file."""
    sample_name = re.sub('\.txt$', '', filename)
    sample_name = re.sub('^PCA of ', '', sample_name)
    d = pandas.read_table(PCA_dir + filename)
    # get indexes of bases in clone
    pca_seq = ''.join(d['Clone sequence'].fillna(' '))
    pos = get_pos_in_clone_seq(clone_seq, pca_seq)
#    positions = num_matching_pos(pca_seq, clone_seq)
    d['Position in clone'] = pos
    # FIXME add some measure of alignment quality?
    num_missing = d['Position in clone'].isna().sum()
    print(f'{sample_name} : {num_missing} sites missing position in clone')
    d['Sample name'] = sample_name
    d['Dist using PC1 and PC2'] = np.sqrt( d.PC1**2 + d.PC2**2 ).round(5)
    return d

per_base_PCA = [parse_PCA(f) for f in os.listdir(PCA_dir)]
per_base_PCA = pandas.concat(per_base_PCA)

per_base_PCA.to_csv('PCA_base_filter.csv.gz')

