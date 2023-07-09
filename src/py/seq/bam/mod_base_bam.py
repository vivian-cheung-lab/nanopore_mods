# Writes out modified bases in BAM format.

import array
import pdb

import numpy as np

import pysam

class ModBaseWriter:
    """Writes a BAM file of bases likely to be modified.

    Uses the MM/ML modified base tags described in
    https://samtools.github.io/hts-specs/SAMtags.pdf
    """

    def __init__(self, canonical_base=None):
        """Constructor.

        canonical_base: the canonical base to set this to
            (currently only supports one canonical base)
        This has some settings which can be modified:
        set_canonical_base (default True): should the read's sequence
            be set to the given base?
        """
        self.set_canonical_base = True
        # for now, always using a fixed "canonical" base
        self.canonical_base = canonical_base
        # ??? eventually use this to pick a "canonical" base
        # self.mod_to_canonical = {'c': 'C', 'f': 'C', 'm': 'C'}

    def get_modified_locations(self, read, mods):
        """Gets which bases are modified.

        read: the read, as an AlignedSegment
        mods: the modifications, as a list of (genomic_pos, mod, prob) tuples
        Returns: a list of modifications which were present in the read,
            but as a list of (query_pos, mod, prob) tuples
        """
        # Get aligned (query pos, reference pos) pairs.
        # (Presumably these are sorted.)
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        # create mapping from reference position (now 1-based)
        # to query position
        ref_to_query = {ref_pos+1: query_pos
            for (query_pos, ref_pos) in aligned_pairs}
        # Find which of those are considered modified, and tack on
        # probabilities for them.
        mods = [(ref_to_query[m[0]], m[1], m[2]) for m in mods
            if m[0] in ref_to_query]
        return mods

    def set_bases_to_canonical(self, read, query_locs):
        """Set modified bases to the canonical base.

        This may be necessary, since some tools (such as IGV) don't
            seem to support 'N' as the canonical base.
        read: the read
        query_locs: locations to set to the canonical base
        Side effects: sets the bases at those indices in the read to
            the canonical base. 
        Returns: offsets of modified bases (only counting the canonical base)
        """
        # save qualities (since modifying the query string invalidates them)
        qualities = read.query_qualities
        s = np.array(list(read.query_sequence))
        s[query_locs] = self.canonical_base
        # fill in (modified) sequence, and qualities
        read.query_sequence = ''.join(s)
        read.query_qualities = qualities

    def get_mod_tag(self, read, mods, mod_name):
        """Gets tag info for one modification.

        read: the read
        mods: list of modifications (as for set_mod_tags)
        mod_name: name of the modification to include
        Returns: a pair, of the MM tag (as a string),
            and the ML tag (as a numpy array)
        """
        # get just these modifications
        mods = [m for m in mods if m[1] == mod_name]
        # get query indices
        query_locs = np.array([m[0] for m in mods])
        # get indices of modified characters
        s = np.array(list(read.query_sequence))
        base_counts = np.cumsum( s == self.canonical_base )
        # get indices (only counting the canonical base)
        i = base_counts[query_locs]
        offsets = np.concatenate([[i[0]-1], np.diff(i)-1])
        offset_string = ','.join([str(i) for i in offsets])
        return (f'{self.canonical_base}+{mod_name},{offset_string};',
            np.array([m[2] for m in mods]))

    def set_mod_tags(self, read, mods):
        """Adds modification info to a read.

        Note that this may update the read sequence in order to set the
            "canonical" base there.
        read: an AlignedSegment object
        mods: the modifications, as a list of (pos, mod, prob) tuples, where
            pos: the genomic position of the modification
            mod: the name of the modification
            prob: the probability of the modification
        Currently this assumes that each position only has one modification.
        Side effects: sets the MM and ML tags.
        """
        # look up indices in the read which were modified
        mods = self.get_modified_locations(read, mods)
        # ??? are these sorted? I think so
        # if no bases were modified, don't add tag
        if not mods:
            return read
        # set bases which are modified to the "canonical" base
        query_locs = np.array([m[0] for m in mods])
        self.set_bases_to_canonical(read, query_locs)
        # get list of all the modifications
        mod_names = list(set([m[1] for m in mods]))
        mod_names.sort()
        # get tag data for each of the modifications
        tags = [self.get_mod_tag(read, mods, mod_name)
            for mod_name in mod_names]
        MM_tag = ''.join([t[0] for t in tags])
        ML_tag = np.concatenate([t[1] for t in tags])
        # fill in tags
        read.set_tag('MM', MM_tag)
        # for ML tag, first clip, then convert to array of uint8
        ML_tag = np.clip(np.trunc(255 * ML_tag), 0., 255.)
        ML_tag = ML_tag.astype('uint8')
        ML_tag = array.array('B', ML_tag)
        read.set_tag('ML', ML_tag)
        return read

