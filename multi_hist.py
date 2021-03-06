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

    def get_ind(self, stat, err_tol=10**-12):
        """
        For a scalar stat, return the index number for its bin. 
        """
        if stat < self.min_val - err_tol or stat > self.max_val + err_tol:
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
        
    def hist_str(self, h_num):
        """
        Turn the kth histogram into a string.
        """
        return ','.join(`v` for v in self.hist[h_num,:])
        
    
    def rates(self, ang, epsilon):
        """
        Get probability mass of samples within epsilon of ang for each stat.
        """
        bin_a = self.get_ind(ang - epsilon)
        bin_b = self.get_ind(ang + epsilon)

        return self.hist[:,bin_a:bin_b+1].sum(axis=1)/float(self.hist[0,:].sum())
