# Aligns nanopore reads.

import argparse  # ??? use this ?
import glob
import gzip
import json
import os
import pdb
import re
import subprocess
import sys

import pandas

import pysam

# index to use for aligning
refIndex = '/home/jtburd/work/seq/align/index/minimap2/hg19_male.mmi'

def getGuppyReads(inputDir, outputFile):
  """Merges reads in a directory into one file.

  inputDir: directory, with guppy output in a subdirectory
  outputFile: name of output .fastq.gz file
  Side effects: writes the merged file .fastq.gz format
  Returns: pathname of the guppy_output/ directory
  """
  guppyOutputDirs = glob.glob(inputDir + '/*/*/guppy_output/')
  # XXX support older name for Guppy output
  if len(guppyOutputDirs) != 1:
    guppyOutputDirs = glob.glob(inputDir + '/*/*/guppy5_output/')
  if len(guppyOutputDirs) != 1:
    print('Couldn\'t find a unique guppy_output/ dir; exiting.')
    sys.exit(1)
  guppyOutputDir = guppyOutputDirs[0]
  passReadsDir = guppyOutputDirs[0] + '/pass/'
  readFiles = glob.glob(passReadsDir + '/*.fastq.gz')
  # only merge files if merged file doesn't exist
  if not os.path.exists(outputFile):
    # just concatenating the gzipped files (this results in
    # somewhat less compression, but is really fast)
    subprocess.run(['cat']+readFiles, stdout=open(outputFile, 'wb'))
  return guppyOutputDir

def alignMinimap2(fastqFile, alnOutputBase):
  """Aligns samples using minimap2.

  fastqFile: the fastq file to read
  alnOutputBase: base name of the output files to write
  Side effects:
  - writes aligned reads to <base>.bam
    If this file exists, returns without doing anything.
  - writes index in <base>.bam.bai
  """
  bamFile = alnOutputBase + '.bam'
  if os.path.exists(bamFile):
    return
  print('[running minimap2]')
  p1 = subprocess.Popen(['minimap2', '-t', '64',
    '-a', '-L',
    # this preset is for RNA
    '-x', 'splice',
    # hopefully the above preset will avoid the need for this
    #  '--secondary=no',
    refIndex, fastqFile],
    stdout=subprocess.PIPE)
  p2 = subprocess.Popen(['samtools', 'sort', '-o', bamFile,
    '-T', 'reads.tmp'], stdin=p1.stdout)
  # ??? I don't completely understand these
  p1.stdout.close()
  p2.communicate()
  subprocess.run(['samtools', 'index', bamFile])

def writeStatistics(guppyOutputDir, alnOutputBase, stats_filename):
  """Writes (very basic!) basecalling and read mapping statistics.

  This makes many assumptions about file layout.
  guppyOutputDir: path to the guppy_output/ dir
  alnOutputBase: base name of output files (including
    <outputBase>.bam)
  stats_filename: file in which to write statistics.
  Side effects: writes a file of stats to statsFile
  """
  # skip this if the output file exists
  if os.path.exists(stats_filename):
    return
  # get read counts
  guppyStats = pandas.read_table(guppyOutputDir + '/sequencing_summary.txt')
  passReads = sum(guppyStats.passes_filtering)
  failReads = sum(1-guppyStats.passes_filtering)
  totalReads = passReads + failReads
  # get mapped read count
  bamFile = alnOutputBase + '.bam'
  samtoolsOutput = subprocess.check_output(['samtools', 'view',
    '-F', str(0x904), '-c', bamFile])
  alignedReads = int(samtoolsOutput)
  # parse the run info file
  run_info_json_file = glob.glob(guppyOutputDir + '/../report*.json')[0]
  # ??? for some reason, this needs to be opened in binary mode
  # I don't know why, but it seems to work if I do, so doing that.
  with open(run_info_json_file, 'rb') as f:
    run_info = json.load(f)
  flow_cell_id = run_info['protocol_run_info']['flow_cell']['flow_cell_id']
  # hopefully this is the number of pores available beforehand
  pores_available = (run_info['acquisitions'][1]['acquisition_run_info']
    ['bream_info']['mux_scan_results'][0]['counts']['single_pore'])
  # this is assumed to be the name for the sample
  sampleName = os.path.basename(alnOutputBase)
  with open(stats_filename, 'w') as statsFile:
    statsFile.write('Sample,Total reads,Passing reads,% passing,Aligned reads,% aligned,Flow cell name,Num pores before\n')
    statsFile.write(','.join([
      sampleName,
      str(totalReads),
      str(passReads), '{:.2%}'.format(passReads / totalReads),
      # ??? should this be a percentage of reads which pass?
      str(alignedReads), '{:.2%}'.format(alignedReads / totalReads),
      flow_cell_id,
      str(pores_available)
    ]) + '\n')

def write_aligned_fastq(bam_file, fastq_in, fastq_out):
    """Filters a .fastq.gz file to those which aligned.

    bam_file: name of the BAM file
    fastq_in: name of the FASTQ file to read
    fastq_out: name of the FASTQ file to write
    Side effects: writes out just the reads which are in the .bam file
    """
    # get set of read names 
    with pysam.AlignmentFile(bam_file) as bam:
        read_names = [r.query_name for r in bam.fetch()]
    read_names = set(read_names)
    # filter reads
    # ??? can the input file be compressed?
    with pysam.FastxFile(fastq_in) as fin, gzip.open(fastq_out, mode='wt') as fout:
        for entry in fin:
            if entry.name in read_names:
                fout.write(str(entry) + '\n')

if __name__ == '__main__':
  cmd = sys.argv[1]
  if cmd == 'minimap2':
    # Map using minimap2.
    # This is assumed to be run in the directory with guppy_output
    # two levels up (that is, './*/*/guppy_output' should match one
    # directory). Output files will be named based on the current
    # working directory.
    outputName = os.path.basename(os.getcwd())
    # the .fastq.gz file will be in this directory
    fastqFile = os.getcwd() + '/' + outputName + '.fastq.gz'
    # concat .fastq.gz files (and find guppy_output/ dir)
    guppyOutputDir = getGuppyReads(os.getcwd(), fastqFile)
    # base name for alignment output
    alnOutputDir = re.sub('\/Fastq\/', '/aln/', os.getcwd())
    os.makedirs(alnOutputDir, exist_ok=True)
    alnOutputBase = alnOutputDir + '/' + outputName
    # align
    alignMinimap2(fastqFile, alnOutputBase)
    # write out various stats
    writeStatistics(guppyOutputDir, alnOutputBase)
    sys.exit(0)

  if cmd == 'writeStats':
    # Write out statistics.
    # This is somewhat deprecated, since it makes many assumptions about
    # output file layout.
    # Again, this is assumed to be run in the directory with guppy_output
    # two levels up (that is, './*/*/guppy_output' should match one
    # directory). Output files will be named based on the current
    # working directory.
    outputName = os.path.basename(os.getcwd())
    # the .fastq.gz file will be in this directory
    fastqFile = os.getcwd() + '/' + outputName + '.fastq.gz'
    # concat .fastq.gz files (and find guppy_output/ dir)
    guppyOutputDir = getGuppyReads(os.getcwd(), fastqFile)
    # base name for alignment output
    alnOutputDir = re.sub('\/Fastq\/', '/aln/', os.getcwd())
    os.makedirs(alnOutputDir, exist_ok=True)
    alnOutputBase = alnOutputDir + '/' + outputName
    # write out stats
    # skip this if the output file exists
    stats_filename = alnOutputBase + '.stats_and_flow_cell.csv'
    # write out statistics
    writeStatistics(guppyOutputDir, alnOutputBase, stats_filename)
    sys.exit(0)

  print('unknown command: ' + cmd)
  sys.exit(1)

