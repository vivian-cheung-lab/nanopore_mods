#!/usr/bin/env python3
# Writes BAM files including modified bases.
# FIXME
# - write out track of which sites have a large "distance"?
#   (or else color modifications by distance)

import os
import pdb
import re
import subprocess
import sys

import numpy
import pandas

import pysam

sys.path.append('../../../../../py/')

import seq.bam.mod_base_bam

# import seq.nanopore.event.event_outlier
import seq.nanopore.squiggle.read_event

data_dir = '../../polish/f5c/IVT/mods/Snakemake/'

# table of sample names
IVT_samples = pandas.read_csv(data_dir + 'IVT_samples.csv')

ref_fasta_file = '../../../../../../chr19_AANCR_clone_only_230411.fa'

def get_read_event_data_by_sample_name(sample_name):
    """Gets read events based on the 'sample name'. """
    return seq.nanopore.squiggle.read_event.ReadEventData(
        data_dir + 'fastq/' + sample_name + '.fastq.gz',
        data_dir + 'events/' + sample_name + '.mergedEvent.tsv.gz')

def get_read_event_data(name, concentration):
    """Gets read events from the name of the modification, and its concentration."""
    s = IVT_samples[ (IVT_samples.short_name==name)
                   & (IVT_samples.mod_concentration==concentration) ]
    # print(s)
    return get_read_event_data_by_sample_name(s.FASTQ_file_base_name.iloc[0])

pca_base_filter = pandas.read_csv('PCA_base_filter.csv.gz')

class ModCaller:
    """Calls modified sites in reads.

    This is based on mismatches and dwell time.
    """

    def __init__(self, sample_name, mod_name, ref_fasta_file, canonical_base, mismatch_bases, dwell_z_cutoff=None):
        """Constructor."""
        self.canonical_base = canonical_base
        self.mismatch_bases = mismatch_bases
        self.dwell_z_cutoff = dwell_z_cutoff
        # connection to the reads
        self.bam = pysam.AlignmentFile(data_dir + 'bam_no_splicing/'
            + sample_name + '.bam')
        # get connection to read data (for getting event stats)
        self.read_events = get_read_event_data_by_sample_name(sample_name)
        # get dwell times for this sample
        self.dwell_times = self.get_dwell_time_table()
        # get control sample dwell time stats
        self.control_event_stats = pandas.read_csv(
            data_dir + 'event_stats/chr19:45406985-45408892/220822_IVT_AANCR_RNA_nobiotinC.csv.gz')
        # get central base
        self.control_event_stats['clone_base'] = (
            self.control_event_stats.reference_kmer.str[2])
        self.control_event_stats.set_index('pos1', inplace=True)
        # the reference sequence
        ref_fasta = pysam.FastaFile(ref_fasta_file)
        clone_seq = ref_fasta.fetch('chr19', 45406985-1, 45408892)
        # indices which were at the canonical site in the clone seq
        # self.site = set(self.control_event_stats.clone_base.loc[self.control_event_stats.clone_base==canonical_base].index.tolist())
        self.site = set([i + 45406985 for i in range(len(clone_seq))
            if clone_seq[i]==canonical_base])
        # get table of modified bases
        self.get_pca_filter(mod_name)
        # number of reads to include
        self.num_reads = 20

    def get_pca_filter(self, mod_name):
        """Gets the PCA filter of bases for this sample."""
        # first, subset to this sample
        sample_name_1 = mod_name + ' 1:2'
        if mod_name == 'biotin-C':
            sample_name_1 = 'biotin 1:2'
        if mod_name == 'Y':
            sample_name_1 = 'pseudoU 1:2'
        if mod_name == 'm1A':
            sample_name_1 = 'm1A'
        self.pca_base = pca_base_filter[ pca_base_filter['Sample name'] == sample_name_1 ]
        # print various relevant stats
        print(f'{sample_name_1}')
        print(f'Number of canonical base {self.canonical_base} = {len(self.site)}')
        self.pca_base = self.pca_base[ ~ self.pca_base['Dist using PC1 and PC2'].isna() ]
        print(f'Number of rows with distance defined = {self.pca_base.shape[0]}')
        # self.pca_base = self.pca_base[ self.pca_base['Dist using PC1 and PC2'] >= 1 ]  # was 2.1497
        # print(f'Number of those beyond cutoff = {self.pca_base.shape[0]}')

        # restrict to sites with the canonical base
        self.pca_base = self.pca_base[ ~ self.pca_base['Position in clone'].isna() ]
        self.pca_base['pos'] = (45406984 + self.pca_base['Position in clone']).astype(int)
        self.pca_base = self.pca_base[ [p in self.site for p in self.pca_base.pos] ]
        print(f'Number of those which have the canonical base = {self.pca_base.shape[0]}')

        # sort by distance (decreasing), and pick furthest away
        self.pca_base = self.pca_base.sort_values('Dist using PC1 and PC2', ascending=False)
        cutoff = 0.33333333333333333   # was 0.2
        self.pca_base = self.pca_base[ :int(cutoff * self.pca_base.shape[0]) ]
        # write out BED file, to show distances
        self.write_PCA_BED_file(mod_name)
        # convert this to a set
        self.pca_filter = set(self.pca_base.pos.tolist())
        print(f'Number of those near cutoff = {len(self.pca_filter)}')

    def write_PCA_BED_file(self, mod_name, output_dir='.'):
        """Writes a BED file of modified sites, based on the PCA."""
        # assemble columns, in BED format
        BED_file = pandas.DataFrame({
            'chrom': 'chr19',
            'start': self.pca_base.pos - 1,
            'end': self.pca_base.pos,
            'name': self.pca_base['Dist using PC1 and PC2'],
            'score': self.pca_base['Dist using PC1 and PC2'],
            'strand': ''})
        # sort, and write out file
        BED_file = BED_file.sort_values('start')
        BED_file.to_csv(output_dir + '/' + mod_name + '.bed', sep='\t',
            header=False, index=False)

    def get_dwell_time_table(self):
        """Gets table of dwell times for this sample."""
        events = self.read_events.get_events_at_region('chr19:45406985-45408892')
        # for practice, a subset of events
        # events = self.read_events.get_events_at_region('chr19:45406985-45407185')
        events = events[['read_name', 'pos1', 'dwell_time']]
        events.set_index(['read_name', 'pos1'], inplace=True)
        return events

    def get_mod_score_from_dwell(self, read_name, ref_pos):
        """Gets 'modified' score based on dwell time."""
        # if there's no z cutoff, assume "not modified"
        if not self.dwell_z_cutoff:
            return None
        try:
            # get background stats
            bg = self.control_event_stats.loc[ref_pos]
            # get dwell time for this event
            dwell_time = self.dwell_times.loc[read_name, ref_pos].dwell_time
            # compute z-score
            z = (dwell_time - bg.dwell_time_mean) / bg.dwell_time_std
            # check if the z is beyond the cutoff (either higher or lower)
            if (self.dwell_z_cutoff > 0) and (z >= self.dwell_z_cutoff):
                return 1.
            if (self.dwell_z_cutoff < 0) and (z <= self.dwell_z_cutoff):
                return 1.
            # if we reach here, the z-score wasn't too extreme
            return None
        except KeyError:
            return None

    def get_modifications_one_read(self, read):
        """Gets modifications for one read."""
        # this will store the modified locations
        mod_scores = []
        # get aligned (query pos, reference pos) pairs
        pairs = read.get_aligned_pairs(matches_only=True)
        # loop through these
        for (query_pos, ref_pos) in pairs:
            # ref_pos is 0-based! convert it to 1-based
            ref_pos += 1
            # only call a modification if it's the relevant base
            if not ref_pos in self.pca_filter:       # was self.site:
                continue
            base_call = read.query_sequence[query_pos]
            is_mismatch = (base_call in self.mismatch_bases)
            is_high_z = self.get_mod_score_from_dwell(read.query_name, ref_pos)
            # possibly call this as a modification; here, we're just using
            # arbitrary modifications to color bases
            # ... except, now trying to color by score
            score = None
#            if is_mismatch and not is_high_z:
#                score = 1./2.
#            if not is_mismatch and is_high_z:
#                score = 1./2.
            if is_mismatch or is_high_z:
                score = 1.
            if score:
                mod_scores += [(ref_pos, 'm', score)]
        return mod_scores

    def get_modifications(self):
        """Gets modifications for the reads."""
        mods = {}
        for read in self.bam.fetch():
            mods[read.query_name] = self.get_modifications_one_read(read)
            if len(mods) >= self.num_reads:
                break
        return mods

def write_BAM_with_mod_calls(mod_name, concentration,
        mod_base, mismatch_bases, dwell_z_cutoff):
    """Writes out a file with modification calls."""
    s = IVT_samples[ (IVT_samples.short_name==mod_name)
            & (IVT_samples.mod_concentration==concentration) ]
    sample_name = s.FASTQ_file_base_name.iloc[0]
    output_bam_file = sample_name + '_with_mods.bam'
    # skip if output file already exists
    # if os.path.exists(output_bam_file):
    #     return
    print(f'{mod_name} {concentration}')
    input_bam_file = data_dir + '/bam_no_splicing/' + sample_name + '.bam'
    # call modifications
    mod_caller = ModCaller(sample_name, mod_name, ref_fasta_file, mod_base, mismatch_bases,
        dwell_z_cutoff=dwell_z_cutoff)
    w = seq.bam.mod_base_bam.ModBaseWriter(mod_base)
    # loop through the reads
    with pysam.AlignmentFile(input_bam_file) as input_bam, \
            pysam.AlignmentFile(output_bam_file, 'wb', template=input_bam) \
            as output_bam:
        read_count = 0
        for read in input_bam.fetch():
            mods = mod_caller.get_modifications_one_read(read)
            w.set_mod_tags(read, mods)
            output_bam.write(read)
            read_count += 1
            if read_count >= mod_caller.num_reads:
                break
    pysam.index(output_bam_file)

# hack, for purposes of coloring bases as modified:
# red:    m ("5mC")
# green:  c ("5caC")
# yellow: f ("5fC")
# also, orange: g ("5hmU")

# write out the BAM files
z_cutoff = 0.5       # 1.96 would be 95%; 1.65 would be 90%

# the biotin-C file is short, and so runs faster...
write_BAM_with_mod_calls('biotin-C', '1:2', 'C', 'AGT', z_cutoff)

write_BAM_with_mod_calls('Y', '1:2', 'T', 'ACG', z_cutoff)

write_BAM_with_mod_calls('m5C', '1:2', 'C', 'AGT', z_cutoff)
write_BAM_with_mod_calls('hm5C', '1:2', 'C', 'AGT', z_cutoff)

write_BAM_with_mod_calls('m1A', '1:2', 'A', 'CGT', -z_cutoff)
write_BAM_with_mod_calls('m6A', '1:2', 'A', 'CGT', z_cutoff)

# put them all in a folder
subprocess.run('(cd ..; zip -r IGV_1.zip IGV)', shell=True)

