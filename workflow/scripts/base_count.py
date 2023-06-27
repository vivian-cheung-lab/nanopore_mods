#!/usr/bin/env python3
# Counts bases, to find RDDs.
# Modified version to only deal with one file.
# usage:
# base_count.py [input_file.bam] [output_file.csv.gz] [reference.fa]

# FIXME
# - extend to work with multiple regions
# - get set of regions from input sequence

import argparse
import collections
import io
import multiprocessing
import os
import pdb
import re
import sys

import numpy as np
import pandas
import pybedtools
import pybedtools.featurefuncs
import pysam

parser = argparse.ArgumentParser(
    description='Counts bases at each location in a BAM file.')
parser.add_argument('input_bam_filename',
    help='Name of BAM file to count bases in')
parser.add_argument('output_csv_filename',
    help='Name of the CSV file to write '
        '(if this ends in .csv.gz, then output will be compressed)')
parser.add_argument('reference_fasta',
    help='Name of the reference FASTA file')
parser.add_argument('regions',
    help='BED file of regions to include')

args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__)) + '/'

# regions = pybedtools.BedTool('chr19 45406985 45408892 AANCR 0 +', from_string=True)
# shorter region for practice
# regions = pybedtools.BedTool('chr19 45407500 45407525 AANCR 0 +', from_string=True)

# regions on the - strand
# regions = pybedtools.BedTool('chr9 126813687 126813856 LHX2 0 -', from_string=True)
# regions = pybedtools.BedTool('chr4 76807526 76807650 PPEF2 0 -', from_string=True)
# regions = pybedtools.BedTool('chr7 128828405 128828505 foo 0 -', from_string=True)


# for debugging: print information where there are insertions?
print_insertion_details=False

# if output file exists, skip
if os.path.exists(args.output_csv_filename):
    print(f'output file {args.output_csv_filename} exists; exiting')
    sys.exit(0)

# various utilities

def BedTool_copy(bedtool):
  """Makes an in-memory BedTool."""
  return pybedtools.BedTool(BedTool_to_list(bedtool))

def BedTool_to_list(bedtool):
  """Converts a BedTool to a list.

  bedtool: a BedTool object
  Returns: that, as a list of intervals
  """
  return [i for i in bedtool]

def remove_chr(r):
  """Utility to remove 'chr' from a region name."""
  r.chrom = re.sub('^chr', '', r.chrom)
  return r

# for complementing DNA bases (and some related characters)
# note: doesn't work for RNA
complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
    'I':'I', 'D':'D', 'i':'i', 'X':'x'}

def get_regions(bam_file):
    """Gets regions at which to count bases.

    """
    # get regions bounded by the reads (ignoring indels and splicing)
    reads = pybedtools.BedTool(bam_file).bam_to_bed()
    # merge the regions (paying attention to strand)
    regions = reads.merge(s=True)
    # the strand seems to be in the name, so set strand based on name
    def f(r):
        r = pybedtools.featurefuncs.extend_fields(r, 6)
        r.strand = r.name
        return r
    regions1 = regions.each(f)
    regions1 = BedTool_copy(regions1)
    return regions1

def get_base_counts_like_IGV(column, stranded=True):
  """Gets the base (and indel) counts at a position

  column: a column, as returned by AlignmentFile.pileup()
  stranded: if True, only add to counts from one strand.
  Returns: a dict of counts, indexed first by strand ('+' or '-'), then
    of counts:
    A, C, G, T, 'I' for insertion, or 'D' for deletion (as in IGV),
    or 'i' for total number of bases inserted,
    or 'X' for total number of deletion starts.
    Note that insertions will be counted at the base _after_ the
    current column, but IGV counts them at the base _before_ the
    current column (in genomic coordinates, regardless of strand).
    This is arguably arbitrary; callers can presumably shift that
    column if necessary.
    FIXME avoid the offset in insertion counts?
  """
  # get collection of PileupRead objects
  pileups = column.pileups
  # initialize the counts (for each strand)
  count = {strand: {b:0 for b in 'ACGTIDiX'} for strand in '+-'}
  # loop through the reads "piled up" here 
  for read in pileups:
    # get the strand
    strand = '-' if read.alignment.is_reverse else '+'
    # if there's a base aligned here, then add to its counts
    i = read.query_position
    if i != None:
      base = read.alignment.query_sequence[i]
      if base in 'ACGT':
        count[strand][base] += 1 
        # possibly add to counts on other strand
        if not stranded:
          count['+' if strand == '-' else '-'][base] += 1 
    # if there's an insertion, add to the counts here
    # (but see above about shifting this)
    if read.indel > 0:
      # count number of times there was an insertion here
      count[strand]['I'] += 1
      # also count the total number of bases inserted
      count[strand]['i'] += read.indel
      if print_insertion_details:
        # N.B. here we print the 1-based coordinate which IGV uses (which
        # differs slightly from the 0-based coordinate in the .bam file)
        print(f'pos={column.reference_pos+2} indel={read.indel} read={read.alignment.query_name}')
    # if the next site is a deletion,
    # add to count of deletions starting here
    if read.indel < 0:
        count[strand]['X'] += 1
    # count if this base is a deletion (but not a "refskip", which
    # includes deletions and introns)
    if read.is_del and not read.is_refskip:
      count[strand]['D'] += 1
  return count

def get_base_counts_at_one_region(bam_filename, region, dna_mode=False):
  """Gets counts of each base in a region.

  This uses pysam. It may therefore be somewhat slower, but hopefully
    is easier to understand.
  bam_filename: name of a .bam file
  region: a pybedtools Interval, giving the interval in question
  dna_mode: enables various tweaks for counting DNA bases
  Returns: a pandas object with the genomic position, and base counts
  FIXME avoid the offset in insertion counts?
  """
  bam = pysam.AlignmentFile(bam_filename)
  # these will store the results for each strand
  count_by_strand = {'+':[], '-':[]}
  # get iterator to the pileup (which depends on the mode)
  if dna_mode:
    pileup = bam.pileup(region.chrom, region.start, region.end,
      min_base_quality=0, truncate=True, add_indels=True, stepper='nofilter')
  else:
    # FIXME find default for "stepper", to simplify this
    pileup = bam.pileup(region.chrom, region.start, region.end,
      min_base_quality=0, truncate=True, add_indels=True)
    # XXX for testing
    # stepper='nofilter')
  # loop through a pileup of the reads at this region
  for column in pileup:
    base_count = get_base_counts_like_IGV(column, stranded=(not dna_mode))
    # base position is 0-based, but write it as 1-based
    count_by_strand['+'].append([column.reference_name, column.reference_pos+1, '+'] +
      [base_count['+'][b] for b in 'ACGTIDiX'])
    # for '-' strand, we complement the base
    count_by_strand['-'].append([column.reference_name, column.reference_pos+1, '-'] +
      [base_count['-'][b] for b in 'TGCAIDiX'])
  # Converts the list of tuples to a DataFrame.
  def tuples_to_DataFrame(r):
    counts = pandas.DataFrame.from_records(r,
      columns=['chrom', 'pos', 'strand',
        # XXX using T consistently for now
        'A', 'C', 'G', ('T' if dna_mode else 'U'),
        'I', 'D', 'i', 'X'])
    # Shift offset of I, to (hopefully) agree with IGV. Note that this
    # assumes that all consecutive positions have been included.
    # presumably the first aligned base is never an insertion?
    def shift_bases(x):
      return [0] + x[:-1].to_list()
    counts['I'] = shift_bases(counts['I'])
    # do similarly for 'i', which is count of total bases inserted
    counts['i'] = shift_bases(counts['i'])
    return counts
  r = pandas.concat(
    [tuples_to_DataFrame(count_by_strand[strand]) for strand in '+-'])
  return r

def get_counts_at_bases(bam_filename, regions, dna_mode=False):
  """Gets base and mismatch counts at some regions.

  This adds a slight margin, in order to get the same insertion
    counts as IGV. (Basically, it's a workaround for the fact that
    IGV counts insertions differently from pysam.)
  bam_filename: name of the file
  regions: the regions at which to get counts
  dna_mode: enables various tweaks for counting DNA bases
  Returns: a pandas data frame of counts
  """
  print('getting counts for ' + bam_filename)
  r = []
  with pysam.AlignmentFile(bam_filename) as bam:
    for region in regions:
      # ??? "progress printing" isn't working
      # print('\r' + str(region), end='')
      # add 1nt margin on each side
      region_with_margin = pybedtools.Interval(
        region.chrom, region.start-1, region.end+1, strand=region.strand)
      r1 = get_base_counts_at_one_region(bam_filename,
        region_with_margin, dna_mode=dna_mode)
      # just keep the counts at that base, on that strand
      # r1 = r1.loc[(r1.chrom==region.chrom)
      #   & (r1.pos==region.end)
      #   & (r1.strand==region.strand)]
      r1 = r1[ r1.strand==region.strand ]
      # if this is non-empty, add to the list
      if not r1.empty:
        r.append(r1)
  # if no regions have any bases, fail with an error message
  if not r:
    print('no bases found; exiting')
    sys.exit(1)
  counts = pandas.concat(r)
  # XXX add 'chr' to chromosome name, if it's absent
  # counts.chrom = ['chr' + s if not s.startswith('chr') else s for s in counts.chrom]
  counts = counts.set_index(['chrom', 'pos', 'strand'])
  return counts

def blankColumn(name):
  """Creates a named, blank column (for labeling tables).
  """
  return pandas.DataFrame({name: [None]})

# regions = get_regions(args.input_bam_filename)
regions = pybedtools.BedTool(args.regions)

regions_without_chr = BedTool_copy(regions.each(remove_chr))

def sequence_as_dataframe(fasta_file, regions, column_name):
    """Gets sequence from a fasta file, as a pandas DataFrame.

    fasta_filename: name of the .fasta file
    regions: a BedTools object, using 1-based numbering
    column_name: column name to use for the sequence
    Returns: a pandas DataFrame, indexed by 'chrom' and 'pos',
        with the base as a column named 'column_name'
    """
    fasta = pysam.FastaFile(fasta_file)
    sequences = []
    for region in regions:
        seq = fasta.fetch(region.chrom, region.start-1, region.end).upper()
        # possibly complement (but don't reverse) sequence
        if region.strand=='-':
            seq = [complement[s] for s in seq]
        # convert to RNA alphabet
        seq = [s.replace('T', 'U') for s in seq]
        r = pandas.DataFrame({'chrom': region.chrom,
            'pos': range(region.start, region.end+1),
            'strand': region.strand,
            column_name: list(seq)})
        sequences.append(r)
    r = pandas.concat(sequences)
    return r.set_index(['chrom', 'pos', 'strand'])

def get_counts(bam_file, sample_name):
    """Gets counts for one file."""
    counts = get_counts_at_bases(bam_file, regions)
    counts = counts.reset_index()
    counts['Sample'] = sample_name 
    counts = counts.set_index(['chrom', 'pos', 'strand'])
    return counts

# get counts for this .bam file
# (including a name based on the filename)
counts = get_counts(args.input_bam_filename, os.path.basename(args.input_bam_filename))

# add up event counts
counts['n'] = counts.loc[:,['A','C','G','U','I','D']].sum(axis=1)
counts = counts[counts['n']>0]
# tack on reference sequence
# (omitting clone sequence for now, as we won't always have that)
ref_seq = sequence_as_dataframe(args.reference_fasta, regions, 'Ref')
counts = counts.join(ref_seq)
counts = counts[ ~ pandas.isna(counts.Ref) ]

# compute mismatches
counts1 = counts.reset_index()
mismatches = [counts1.iloc[i].n - counts1.loc[i, [counts1.iloc[i].Ref]][0] for i in range(len(counts1))]
counts['Mismatches'] = mismatches

# write out all the counts
counts.to_csv(args.output_csv_filename)

# add up totals, for each base in the reference sequence
# total_counts = counts.groupby('Ref')[['A','C','G','U','I','D','n']].sum()

