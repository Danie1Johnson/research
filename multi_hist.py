import numpy as np


class MultiHist:
    """
    Array of histograms dynamically updated.
    """

    def __init__(self, min_val, max_val, num_stats, num_bins, strict=True):
        self.min_val = min_val
        self.max_val = max_val
        self.num_stats = num_stats
        self.num_bins = num_bins
        self.hist = np.zeros((num_stats, num_bins + 2))
        #self.unders = np.zeros(num_stats)
        #self.overs = np.zeros(num_stats)
        #self.strict = strict
        self.midpoints = (min_val[stat_num] 
                          + (0.5 + np.arange(-1, self.num_bins + 1))
                          *(self.max_val - self.min_val)/self.num_bins)

    def get_inds(self, stats):
        """
        Take ndarray of length num_stats and return the bin index for each in a ndarray of the same shape.
        """
        inds = np.zeros(self.num_stats)
        for k in range(self.num_stats):
            inds[k] = self.get_ind(stats[k], k)
        return inds

    #def get_ind(self, stat):
    #    """
    #    For a scalar stat, return the index number for its bin. 
    #    """
    #    if stat < self.min_val or stat > self.max_val:
    #        raise Exception("ERROR: stat outside of histogram range " 
    #                        + str(self.min_val) + " <  " + str(stat) + " < " + str(self.max_val))
    #    return int(self.num_bins*(stat - self.min_val)/(self.max_val - self.min_val))


    def get_ind(self, stat):
        """
        For a scalar stat, return the index number for its bin. 
        """
        if stat < self.min_val:
            return 0
        if stat > self.max_val:
            return self.num_bins + 1
        return 1 + int(self.num_bins*(stat - self.min_val)/(self.max_val - self.min_val))

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

    def rates_2(self, nu, stat_num):
        """
        Get probability mass of samples within epsilon of ang for each stat.
        """
        bin_a = 1
        bin_b = self.get_ind(nu)

        #return self.hist[:,bin_a:bin_b+1].sum(axis=1)/float(self.hist[0,:].sum())
        return self.hist[stat_num,bin_a:bin_b+1].sum(axis=1)/float(self.hist[stat_num,:].sum())
