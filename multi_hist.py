import numpy as np


class MultiHist:
    """
    Array of histograms dynamically updated.
    """

    def __init__(self, min_val, max_val, num_stats, num_bins):
        self.min_val = min_val
        self.max_val = max_val
        self.num_stats = num_stats
        self.num_bins = num_bins
        self.hist = np.zeros((num_stats, num_bins))
        self.midpoints = min_val + (0.5 + np.arange(num_bins))*(max_val-min_val)/num_bins 

    def get_inds(self, stats):
        """
        Take ndarray of length num_stats and return the bin index for each in a ndarray of the same shape.
        """
        inds = np.zeros(self.num_stats)
        for k in range(self.num_stats):
            inds[k] = self.get_ind(stats[k])
        return inds

    def get_ind(self, stat):
        """
        For a scalar stat, return the index number for its bin. 
        """
        if stat < self.min_val or stat > self.max_val:
            raise Exception("ERROR: stat outside of histogram range " 
                            + str(self.min_val) + " <  " + str(stat) + " < " + str(self.max_val))
        return int(self.num_bins*(stat - self.min_val)/(self.max_val - self.min_val))

    def add_stats(self, stats):
        """
        Take new stats and update the histogram
        """
        new_inds = self.get_inds(stats)
        for k in range(self.num_stats):
            self.hist[k,new_inds[k]] += 1
        
