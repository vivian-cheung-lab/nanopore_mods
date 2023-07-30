# Computes how much different events are "outliers",
# based on event statistics.

import pdb
import sys

import numpy as np
import pandas

# for Mahalanobois distance
import scipy.spatial.distance

import seq.nanopore.squiggle.read_event

class EventOutlier:
    """Finds bases with events which are outliers."""

    def __init__(self, region_string):
        """Constructor.

        region_string: the region to include
        """
        self.region_string = region_string
        self.region = seq.nanopore.squiggle.read_event.parse_region_string(
            region_string)

    def get_event_stats_as_matrix(self, read_data):
        """Gets event stats at a region.

        This assumes the region is on the + strand.
        read_data: object for getting the read data
        Returns a dict with keys:
            read_names: the names of the reads
            event_stats: a Numpy array with indices
                read_number: which read this was
                base: the index of the base in the region
                event_type: 0 for current (pA), 1 for dwell time (seconds)
        """
        # get events here
        events = read_data.get_events_at_region(self.region_string)
        # reshape into a table
        events_table = events.pivot(columns='pos1', index='read_name',
            values=['event_level_mean', 'dwell_time'])
        # convert into a numpy object
        region_length = self.region['end'] - self.region['start']
        x = np.full((events_table.shape[0], 2, region_length), np.nan)
        j = events_table.loc[:,'event_level_mean'].columns.to_numpy() - self.region['start']
        x[:, 0, j] = events_table.loc[:, 'event_level_mean']
        x[:, 1, j] = events_table.loc[:, 'dwell_time']
        return {'read_names': events_table.index.to_list(),
                'event_stats': x}

    def make_normal_model(self, read_data):
        """Makes a Gaussian model of event stats.

        Assume num_sites is the number of sites.
        Returns a dict, with keys:
            n, of shape (num_sites): the number of observations for this site
            mu, of shape (num_sites, 10): the mean at this site
            p, of shape (num_sites, 10, 10): the precision at this site
            (When one of these sites is missing, this will have n=0,
                and NaN for all of mu and p.)
        """
        num_sites = self.region['end'] - self.region['start']
        r = {'n': np.zeros(num_sites),
            'mu': np.full((num_sites, 10), np.nan),
            'p': np.full((num_sites, 10, 10), np.nan)}
        event_stats = read_data['event_stats']
        # ??? log-transform the dwell times?
        # these numbers are for events centered at these locations, so we
        # omit the first two sites, and last two sites
        for j in range(2, num_sites-2):
            # get events in this window, as a matrix with one row per read,
            # and ten columns (five of event current, five of dwell time)
            x = event_stats[:,(i-2):(i+3),:].reshape((-1, 10))
            # find which rows which have all non-missing data
            i = ~np.isnan(x).all(axis=1)

            # ...


        pass


    def get_outlier_scores(self, read_data):
        """Gets scores of which reads have events 'further' from the model.



        """
        pass

