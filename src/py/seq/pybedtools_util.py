# Utilities for use with pybedtools.

import pdb
import re

import pybedtools


def BedTool_to_list(bedtool):
  """Converts a BedTool to a list.

  bedtool: a BedTool object
  Returns: that, as a list of intervals
  """
  return [i for i in bedtool]

def BedTool_copy(bedtool):
  """Makes an in-memory BedTool.

  This avoids having a BedTool be "consumed" when reading
  it (which means that it can't be read more than once).
  (It's quite possible that I'm doing something wrong, and
  that this isn't be necessary.)
  """
  return pybedtools.BedTool(BedTool_to_list(bedtool))

def rename_chrom(name_type):
  """Renames 'chrom' field of an Interval.

  name_type: currently either 'chr' or 'no_chr'
  Returns: function which renames the 'chrom' field to
    match name_type.
  For use with pybedtools.each.
  """
  # functions to tweak chromosome names 
  def add_chr(f):
    if not f.chrom.startswith('chr'):
      f.chrom = 'chr' + f.chrom
    return f
  def remove_chr(f):
    if f.chrom.startswith('chr'):
      f.chrom = re.sub('^chr', '', f.chrom)
    return f
  # return the appropriate function
  if name_type == 'chr':
    return add_chr;
  if name_type == 'no_chr':
    return remove_chr;
  # FIXME throw an exception if name_type isn't known

