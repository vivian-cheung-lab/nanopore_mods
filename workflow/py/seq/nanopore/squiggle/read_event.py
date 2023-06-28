# Gets read event and raw "squiggle" data.
# FIXME
# - in get_events_at_region(), should filter by strand

import io
import pdb
import re

import pandas
import pysam

from ont_fast5_api.fast5_interface import get_fast5_file

def parse_region_string(region):
    """Parses a region string (as used by UCSC, samtools, rtracklayer, etc.)

    region: a region string, in the format "chrom:start-end" or
        "chrom:start-end:strand"
    Returns: a dict with keys 'chrom', 'start', 'end', and 'strand',
        or None on failure
    """
    # parse including strand
    m = re.match('(.+):(\d+)-(\d+):([+-])$', region)
    if m:
        return {'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
            'strand': m.group(4) };
    # parse not including strand
    m = re.match('(.+):(\d+)-(\d+)$', region)
    if m:
        return {'chrom': m.group(1), 'start': int(m.group(2)), 'end': int(m.group(3)),
            'strand': None };
    return None

class ReadEventData:
    """Gets read event data (at both event and squiggle levels).

    """

    def __init__(self, fastq_file, event_file):
        """Constructor.

        This uses the events, as merged by f5c.

        fastq_file: name of the .fastq.gz file of reads, which should be
            indexed by f5c. (The file mapping read names to .fast5
            filenames is assumed to be named this + '.index.readdb' .)
        event_file: name of the file of merged events from f5c.
            This should be tabix-compressed and indexed.
        """
        self.fastq_file = fastq_file
        self.event_file = event_file
        self.read_name_to_fast5_file = fastq_file + '.index.readdb'
        # we read the mapping from read name to fast5 file on-demand;
        # since it's sorted by read name, we could also use binary search
        self.read_name_to_fast5 = None
        # open the event file
        self.events_tabix = pysam.TabixFile(event_file)
        # sample rate XXX for now, we assume this is constant.
        # FIXME when constructor is called, open an arbitrary
        # fast5 file, and get this information from there?
        self.sample_rate = 3012.

    def load_read_name_to_fast5(self):
        """Reads in the mapping from read name to fast5 file name."""
        # if this is already read in, skip this
        if self.read_name_to_fast5 is not None:
            return
        self.read_name_to_fast5 = pandas.read_table(
            self.read_name_to_fast5_file, header=None)
        self.read_name_to_fast5.columns = ['read_name', 'fast5_file']
        self.read_name_to_fast5.set_index('read_name', inplace=True)

    def get_events_at_region(self, region, include_dwell_time=True):
        """Gets the events at a region.

        region: the region, as a Samtools-style region string
        include_dwell_time: if True, include the "dwell time"
            (in seconds) as an additional column
        Returns: a pandas DataFrame of the events in that region
        """
        region1 = parse_region_string(region)
        if not region1:
            return None
        # we fetch sites with genomic position shifted by 3nt,
        # so that the event centers are in the region
        events1 = list(self.events_tabix.fetch(region1['chrom'],
            max(region1['start']-3, 0), region1['end']-3))
        # convert list of rows (each of which is tab-separated)
        # to one long string
        events1 = '\n'.join(events1) + '\n'
        # parse this string
        events = pandas.read_table(io.StringIO(events1),
            names = ['contig', 'position', 'reference_kmer',
                'read_name', 'strand', 'event_index',
                'event_level_mean', 'event_stdv', 'event_length',
                'model_kmer', 'model_mean', 'model_stdv',
                'standardized_level', 'start_idx', 'end_idx'])
        # the position of the central base of the listed 5-mer
        events['pos1'] = events.position + 3
        # the strand
        events['strand1'] = ['+' if eq else '-'
            for eq in (events.model_kmer==events.reference_kmer)]

        # filter out reads with multiple events at a given position
        events = events.drop_duplicates(
            subset=['contig', 'pos1', 'strand', 'read_name'], keep=False)

        # FIXME filter by strand
        # the central base (in model_kmer; so taking strand into account)
        events['central_base'] = events.model_kmer.str[2]
        # possibly include dwell time
        if include_dwell_time:
            events['dwell_time'] = (
                (events.end_idx - events.start_idx) / self.sample_rate)
        return events

    def get_raw_trace(self, read_name, start, end):
        """Gets the raw trace for a given read_id.

        read_name: name of the read
        start, end: indices in the trace to get
        Returns: the trace data for that read, from start to end
            (or None, if the read couldn't be found)
        """
        # if this mapping isn't present, load it
        if not self.read_name_to_fast5:
            self.read_name_to_fast5 = pandas.read_table(
                self.read_name_to_fast5_file)
            self.read_name_to_fast5.columns = ['read_name', 'fast5_file']
            self.read_name_to_fast5.set_index('read_name', keep=True)
        # get name of the file to read
        try:
            traceFile = self.readIdToFilename.loc[read_id][0]
        except KeyError:
            print('couldn\'t get trace for read_id=' + read_id)
            return None
        # get the traces 
        with get_fast5_file(traceFile, mode='r') as f5:
            read = f5.get_read(read_id)
            info = read.get_channel_info()
            # get the scaled current data
            y1 = read.get_raw_data(start, end, scale=True)
            return np.array(y1)

    def get_raw_traces_at_region(self, region, max_reads=10000000000):
        """Gets a generator for reads at a region.

        Note: this now calls a generator, which isn't as well tested.
        region: the region at which to get events
        max_reads: maximum number of reads to return
        Returns: a list with up to max_reads reads
            (in the format returned by get_raw_traces_at_region_gen)
        """
        # based on
        # https://stackoverflow.com/questions/5234090/how-to-take-the-first-n-items-from-a-generator-or-list
        gen = self.get_raw_traces_at_region_gen(region)
        r = [x for _, x in zip(range(max_reads), gen)]
        return r

    def get_raw_traces_at_region_gen(self, region):
        """Gets all traces at a region (as a generator).

        region: the region at which to get events (based on where
            the center of the event is)
        max_reads: the maximum number of traces to return
        Returns a list of dicts, each with keys:
            read_name: the name of the read
            events: the events for that read (as a pandas DataFrame)
            y: the raw data, as a numpy array
            start, end: the indices included in y
        """
        # make sure this mapping is present
        self.load_read_name_to_fast5()
        # get events to include
        events = self.get_events_at_region(region)
        # for each read, get the name of the fast5 file, and
        # the indices in the raw array to include
        events_by_read_name = events.groupby(['read_name'])
        # FIXME sort this by read_name, and limit number of reads here,
        # so that it always returns the reads for a given region?
        read_info = pandas.DataFrame({
            'start': events_by_read_name['start_idx'].min(),
            'end': events_by_read_name['end_idx'].max()
        })
        # tack on which fast5 file contains each read
        read_info = read_info.join(self.read_name_to_fast5)
        # switch to indexing by fast5_file
        read_info = read_info.reset_index().set_index(['fast5_file'])
        # loop through the fast5 files
        for fast5_file in set(read_info.index):
            # get all the reads in that file
            reads_in_file = read_info[read_info.index==fast5_file]
            # open that file
            with get_fast5_file(fast5_file, mode='r') as f5:
                for (f, read_name, start, end) in reads_in_file.itertuples():
                    read = f5.get_read(read_name)
                    info = read.get_channel_info()
                    # get the scaled current data for that part of the read
                    y = read.get_raw_data(start=start, end=end, scale=True)
                    yield {
                        'read_name': read_name,
                        # also including the events, as they'll probably
                        # usually be needed
                        'events': events[ events.read_name == read_name ],
                        'start': start,
                        'end': end,
                        'y': y
                    }

def summarize_event_file(input_file, output_file, reads_per_sample=100):
    """Writes a summary of the events in a file.

    N.B.: this is very specific to the AANCR sequence!
    Hopefully this will be convenient for plotting. This:
    - adds on modification and concentration.
    - adds on "position in clone".
    - only includes a fixed number of (non-missing) events per sample.
    input_file: the events file to read
    output_file: the (sampled) events file to write
    num_reads: how many reads to include
    Returns: a pandas DataFrame of the file's contents
    """
    events = pandas.read_table(input_file,
        names = ['contig', 'position', 'reference_kmer',
            'read_name', 'strand', 'event_index',
            'event_level_mean', 'event_stdv', 'event_length',
            'model_kmer', 'model_mean', 'model_stdv',
            'standardized_level', 'start_idx', 'end_idx'])
    # filter out missing events
    events = events[ events.model_kmer != 'NNNNN' ]
    # count events per read
    events_per_read = events.groupby('read_name', as_index=False).size()
    events_per_read = events_per_read.sort_values(
        ['size', 'read_name'], ascending=[False, True])
    # select reads with the most events
    reads_with_most_events = events_per_read[:reads_per_sample].read_name
    events = events.merge(reads_with_most_events)
    # add on various annotation
    # the position of the central base of the listed 5-mer
    events['pos1'] = events.position + 3
    # the strand
    events['strand1'] = ['+' if eq else '-'
        for eq in (events.model_kmer==events.reference_kmer)]
    # convert T's in sequence to U's
    events['model_kmer'] = events.model_kmer.str.replace('T', 'U')
    # the central base (in model_kmer; so taking strand into account)
    events['central_base'] = events.model_kmer.str[2]
    # the position of the base in the clone
    events['pos_in_clone'] = events.pos1 - 45406984
    # XXX round several of these
    events = events.round({
        'event_level_mean': 3,
        'event_stdv': 4,
        'event_length': 6})
    # subset the columns
    events = events[['contig', 'pos1', 'pos_in_clone', 'strand1',
        'model_kmer', 'central_base',
        'event_level_mean', 'event_stdv', 'event_length']]
    # write out the events 
    events.to_csv(output_file, sep='\t')

