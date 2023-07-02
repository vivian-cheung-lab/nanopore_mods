# Polishes nanopore reads.
# Currently uses f5c.
# FIXME:
# - index events using tabix?
# - only store merged events? (this is probably the only info we need)

import glob
import gzip
import os
import pdb
import re
import subprocess
import sys

import pybedtools
import pysam

class Polisher:
  """Polishes nanopore reads using f5c.

  """

  def __init__(self, fast5Dir, mergedFastqGzFile, bamFile,
        eventOutputBase=None):
    """Constructor.

    fast5Dir: name of directory containing .fast5 files
    mergedFastqGzFile: name of .fastq.gz file
    bamFile: file of aligned reads
    eventOutBase: where to write events; this defaults to
        a file in the directory containing the .bam files
    """
    self.fast5Dir = fast5Dir
    self.mergedFastqGzFile = mergedFastqGzFile
    self.bamFile = bamFile
    # basename of where to write events
    if eventOutputBase:
      self.eventOutputBase = eventOutputBase
    else:
      self.eventOutputBase = re.sub('\.bam', '', self.bamFile)
    # number of threads to use, for multithreaded steps
    self.numThreads = 4

  def indexReads(self):
    """Indexes reads using f5c.

    f5c needs to be in $PATH.
    """
    # if output file seems to exist, skip this
    if os.path.exists(self.mergedFastqGzFile + '.index'):
      return
    subprocess.run(['f5c', 'index',
      # f5c recommends not using -s (summary file) option
      # increase number of threads
      '-t', str(self.numThreads), '--iop', str(self.numThreads),
      '-d', self.fast5Dir,
      # the index files will be written in this directory
      # (which seems like a fairly convenient location)
      self.mergedFastqGzFile])

  def write_merged_events(self, genome_fasta):
    """Aligns and merges events.

    This uses f5c's event-merging, and uses bgzip and tabix to
    compress and index the output.
    genome_fasta: name of the reference genome
    Side effects: writes aligned events.
    """
    # the output filename
    event_output_filename = (self.eventOutputBase
      + '.mergedEvent.tsv.gz')
    # if output exists, skip this
    if os.path.exists(event_output_filename):
      return
    print('[running f5c eventalign]')
    p1 = subprocess.Popen(['f5c', 'eventalign',
      '-r', self.mergedFastqGzFile,
      '-b', self.bamFile,
      '-g', genome_fasta,
      '--rna',
      # write coordinates in signal
      '--signal-index',
      # collapse adjacent events
      '--collapse-events',
      # Print read names, since the "read numbers" are confusing.
      # This makes the event file huge; at least it's compressed.
      '--print-read-names',
      '-t', str(self.numThreads),
      # XXX adding region, for practice
      '--iop', str(self.numThreads)],
      stdout=subprocess.PIPE)
    # remove header line
    p2 = subprocess.Popen(['tail', '--lines=+2'],
      stdin=p1.stdout, stdout=subprocess.PIPE)
    # sort events (this seems to already be partly sorted)
    p3 = subprocess.Popen(['sort',
      # N.B. sorting by contig, pos, and read name, in order
      # to remove duplicate lines
      '-k1,1', '-k2,2n', '-k4,4',
      '--parallel=' + str(self.numThreads),
      '--buffer-size=50G'],
      stdin=p2.stdout, stdout=subprocess.PIPE)
    # XXX There are (very occasionally) duplicate rows; removing those.
    p4 = subprocess.Popen(['uniq'],
      stdin=p3.stdout, stdout=subprocess.PIPE)
    # compress events
    p5 = subprocess.Popen(['bgzip',
      # using faster (but less effective) compression
      '--compress-level', '3',
      # ??? not clear how many threads to use
      '--threads', str(self.numThreads)],
      stdin=p4.stdout,
      stdout=open(event_output_filename, 'wb'))
    # wait until these finish
    p1.stdout.close()
    p2.stdout.close()
    p3.stdout.close()
    p4.stdout.close()
    p5.communicate()
    # index output file
    subprocess.run(['tabix', '--zero-based',
      # skip header
      '--skip-lines', '1',
      '--sequence', '1', '--begin', '2', '--end', '2',
      event_output_filename])

  def eventAlign(self, genomeFasta, region=None, eventOutputBase=None):
    """Aligns events.

    genomeFasta: name of the reference genome
    region: region to align (in the format 'chr:start-end')
    eventOutputBase: base name of output (if not given, this
      will be based on .bam file name)
    Side effects: writes aligned events.
    """
    # default to writing output in a file based off of .bam file name
    if not eventOutputBase:
      self.eventOutputBase = re.sub('\.bam', '', self.bamFile)
    else:
      self.eventOutputBase = eventOutputBase
    # name of the output file (possibly including region)
    if region:
      eventOutputFilename = (self.eventOutputBase 
        + '.' + region + '.events.tsv')
    else:
      eventOutputFilename = (self.eventOutputBase + '.events.tsv')
    # if output exists, skip this
    if os.path.exists(eventOutputFilename + '.gz'):
      return
    print('[running f5c eventalign]')
    # optional region option
    regionOption = ['-w', region] if region else []
    subprocess.run(['f5c', 'eventalign',
      '-r', self.mergedFastqGzFile,
      '-b', self.bamFile,
      '-g', genomeFasta,
      '--rna',
      # write coordinates in signal
      '--signal-index',
      # Print read names, since the "read numbers" are confusing.
      # This also makes the event file huge; at least it's gzipped.
      '--print-read-names',
      '-t', str(self.numThreads),
      '--iop', str(self.numThreads)] + regionOption,
      stdout=open(eventOutputFilename, 'wb'))
    # compress events
    # FIXME do this as the file is written?
    # (previous attempts at this resulted in a garbled file)
    # also using lower compression level (but faster, hopefully)...
    subprocess.run(['gzip', '-3', eventOutputFilename])

def filter_events_by_region(region, event_input_file, event_output_file):
  """Filters events to just include those in some regions.

  (This is a temporary hack; hopefully there are better ways
    to get the index.)
  region: BedTool object, giving the regions in question
  event_input_file: input file of events
  event_output_file: output file of events
  Side effects: writes events in those regions to
    event_output_file. (If that file exists, does nothing.)
    (Note that this uses the "pos" field as the event's position.)
  """
  if os.path.exists(event_output_file):
    return
  # make a set containing individual positions to include
  # (this may be expensive in terms of memory, but hopefully is fast)
  region_pos = set()
  for r in region:
    chrom = r.chrom
    for i in range(r.start, r.end):
      region_pos.add(r.chrom + ':' + str(i))
  with gzip.open(event_input_file, 'rt') as event_in, \
      gzip.open(event_output_file, 'wt') as event_out:
# for debugging, don't compress output?
#      open(event_output_file, 'wt') as event_out:
    first_line = event_in.readline()
    event_out.write(first_line)
    # quick check that the file is in the right format
    if not first_line.startswith(
        'contig\tposition\treference_kmer\tread_name\tstrand\tevent_index\t'):
      print('filter_events_by_region: input event file looks wrong')
      sys.exit(1)
    line_number = 0
    for line in event_in:
      a = line.split('\t')
      # ??? should this include strand?
      event_loc = pybedtools.Interval(a[0], int(a[1]), int(a[1]))
      if (a[0] + ':' + a[1]) in region_pos:
        event_out.write(line)
      # XXX these are for debugging
      if line_number and line_number % 100000 == 0:
        print(line_number, end=' ')
        sys.stdout.flush()
        event_out.flush()
        # break    # XXX for debugging
      line_number += 1
  return

def write_bam_matching_reference(ref_fasta_file, bam_in_file, bam_out_file):
    """Writes out reads, with the reference sequence at their location.

    For purposes of training the model, we may need to assume we know what
    the sequence is. Hopefully, this will reduce the number of mismatches
    and indels, and reduce the amount of missing data in polishers' outputs.

    ref_fasta_file: the FASTA reference sequence
    bam_in_file: BAM file to read
    bam_out_file: BAM file to write
    Side effects: writes bam_out_file
    """
    with \
        pysam.FastaFile(ref_fasta_file) as ref_fasta, \
        pysam.AlignmentFile(bam_in_file) as bam_in, \
        pysam.AlignmentFile(bam_out_file, 'wb', template=bam_in) as bam_out:
        # loop through reads in the input file
        for read in bam_in:
            # create a new read, mostly based on the original alignment
            read1 = pysam.AlignedSegment()
            read1.query_name = read.query_name
            read1.flag = read.flag
            read1.reference_id = read.reference_id
            read1.reference_start = read.reference_start
            read1.mapping_quality = read.mapping_quality
            # skip if reference_name isn't present, as presumably
            # the read didn't align
            if not read.reference_name:
                continue
            # get reference sequence at bounds of where that read aligned
            read1.query_sequence = ref_fasta.fetch(
                read.reference_name,
                read.reference_start,
                read.reference_end).upper()
            # set CIGAR string to "all matching"
            read1.cigar = [(pysam.CMATCH, read.reference_length)]
            # pdb.set_trace()
            # this will clear the base qualities, which should be
            # okay in this case 
            bam_out.write(read1)
    pysam.index(bam_out_file)

