import numpy as np


class MultiHist:
    """
    Array of histograms dynamically updated.
    """

    def __inti__(self, min_val, max_val, num_stats, num_bins):
        self.hist = np.zeros((num_stats, num_buns))

    def get_inds(self, stats):
        """
        Take ndarray of length num_stats and return the bin index for each in a ndarray of the same shape.
        """
        inds = np.zeros(self.num_stats)
        for k in range(self.num_stats):
            inds[k] = self.get_ind(stat)
        return inds

    def get_ind(self, stat):
        """
        For a scalar stat, return the index number for its bin. 
        """
        if stat < min_val or stat > max_val:
            raise Exception("ERROR: stat outside of histogram range")
        return int(self.num_bins*(stat - self.min_val)/(self.max_val - self.min_val))

    def add_stats(self, stats):
        """
        Take new stats and update the histogram
        """
        new_inds = self.get_inds(stats)
        for k in range(self.num_stats):
            self.hist[k,new_ind[k]] += 1
        
