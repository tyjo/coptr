"""
coptr_ref.py
======================
Estimate peak-to-trough ratios using complete reference genomes.
"""
import math
import multiprocessing as mp
import numpy as np
import scipy.optimize
import scipy.stats
import sys

from src.print import print_info, print_warning


class QCResult:
    """Data structure to store results of quality checking.

    Parameters
    ----------
        frac_nonzero : float
            Fraction of nonzero bins in 10Kb windows
        nreads : int
            Number of reads after filtering
        passed_qc : bool
            Flag indicating if the sample passed quality thresholds
    """

    def __init__(self, frac_nonzero, nreads, passed_qc):
        self.frac_nonzero = frac_nonzero
        self.nreads = nreads
        self.passed_qc = passed_qc


class ReadFilterRef:
    """Read filtering steps for CoPTR Ref.

    Parameters
    ----------
        min_reads : float
            Minumum read count after filtering
        frac_nonzero : float
            Fraction of nonzero 10Kb bins required
    """

    def __init__(self, min_reads, min_cov):
        self.min_reads = min_reads
        self.min_cov = min_cov


    def filter_reads(self, read_positions, genome_length):
        """Filter out reads in regions of the genome where the coverage
        is either too high or too low.

        Parameters
        ----------
            read_positions : np.array or list of int
                The starting coordinate of each read along the genome
            genome_length : float
                The length of the reference genome

        Returns
        -------
            filtered_read_positions : np.array of int
                The starting position of each read passing filtering criteria
            filtered_genome_length : float
                The length of the genome with filtered regions removed
            passed_qc : bool
                A flag indicating if the sample passed quality control metrics
        """
        read_positions = np.copy(read_positions)
        read_positions, genome_length = self.filter_reads_phase1(read_positions, genome_length)
        filtered_read_positions, filtered_genome_length = self.filter_reads_phase2(read_positions, genome_length)
        qc_result = self.quality_check(filtered_read_positions, filtered_genome_length, self.min_reads, self.min_cov)
        return filtered_read_positions, filtered_genome_length, qc_result


    def quality_check(self, filtered_read_positions, filtered_genome_length, min_reads, frac_nonzero):
        """A basic quality check for required coverage of a genome in a sample.

        Parameters
        ----------
            filtered_read_positions : np.array or list of int
                The starting coordinate of each read along the genome
                after filtering
            filtered_genome_length : float
                The length of the reference genome after filtering
            min_reads : float
                Minumum read count after filtering
            frac_nonzero : float
                Fraction of nonzero 10Kb bins required

        Returns
        -------
            qc_result : QCResult
                Object that stores the result of the quality check
        """
        binned_counts = self.bin_reads_10Kb(filtered_read_positions, filtered_genome_length)
        nonzero_bins = (binned_counts > 0).sum() / binned_counts.size
        nreads = binned_counts.sum()

        passed_qc = True if nonzero_bins >= frac_nonzero and nreads >= min_reads else False

        return QCResult(nonzero_bins, nreads, passed_qc)


    def compute_genomewide_bounds(self, median_read_count):
        """Compute bounds on read counts in 10Kb bins.

        Parameters
        ----------
            median_read_count : float
                The median read count in 10Kb bins

        Returns
        -------
            lower_bound : float
                A lower bound on read counts in 10Kb bins
            upper_bound : float
                An upper bound on read counts in 10Kb bins

        """
        pois_upper = scipy.stats.poisson(2*median_read_count)
        pois_lower = scipy.stats.poisson((1./2)*median_read_count)
        for i in range(2, 17):
            m = 0.5*i
            lower_bound = (1./m)*median_read_count
            upper_bound = m*median_read_count
            if pois_upper.cdf(upper_bound) > 0.99 and (1 - pois_lower.cdf(lower_bound)) > 0.99:
                return lower_bound, upper_bound
        return lower_bound, upper_bound


    def compute_bin_bounds(self, median_read_count):
        """Compute bounds on read counts in a smaller window along the genome.

        Parameters
        ----------
            median_read_count : float
                The median read count in 10Kb bins

        Returns
        -------
            lower_bound : float
                A lower bound on read counts in 10Kb bins
            upper_bound : float
                An upper bound on read counts in 10Kb bins   
        """
        pois = scipy.stats.poisson(median_read_count)
        for i in range(2, 17):
            m = 0.5*i
            lower_bound = (1./m)*median_read_count
            upper_bound = m*median_read_count
            if pois.cdf(upper_bound) - pois.cdf(lower_bound) > 0.99:
                return lower_bound, upper_bound
        return lower_bound, upper_bound


    def compute_rolling_sum(self, read_positions, genome_length):
        """Compute a rolling sum of read counts in 10Kb bins, sliding 1Kb at
        a time.

        Parameters
        ----------
            read_positions : np.array or list of int
                The coordinate of the start positions of each read along the
                genome.
            genome_length : float
                The length of the reference genome

        Returns
        -------
            rolling_counts : np.array of float
                The read count in each rolling bin
            endpoints : np.array of tuple
                Each tuple gives the left (inclusive) and right (exclusive) 
                end point of each bin.
        """
        bin_size = 10000
        step = 1000
        s = read_positions.size

        sorted_reads = np.sort(read_positions)
        endpoints = np.array([(i, i+bin_size) for i in range(0, genome_length, step)])
        rolling_counts = \
            np.array([
                np.searchsorted(sorted_reads, right, side="right") 
                - np.searchsorted(sorted_reads, left, side="right")
                for left,right in endpoints
            ])

        rolling_counts = np.array(rolling_counts)
        endpoints = np.array(endpoints)
        return rolling_counts, endpoints


    def remove_reads_by_region(self, read_positions, genome_length, regions):
        """Remove reads that overlap a region in regions.
        
        Parameters
        ----------
            read_positions : np.array of int
                The coordinate of the start position of each read
            genome_length : int
                The length of the reference genome
            regions : list of tuple of int
                Each tuple gives an interval in the genome to remove

        Returns
        -------
            read_positions : np.array of int
                Update coordiantes for each read after removing the
                regions in region
            new_genome_length : float
                The length of the genome once each region is regions
                is removed.
        """
        new_genome_length = genome_length
        filter_count = 0
        adjustment = 0
        for left,right in regions:
            remove_size = right - left
            new_genome_length -= remove_size

            keep = np.logical_or(read_positions < (left - adjustment), read_positions >= (right - adjustment))
            adjustment += remove_size
            filter_count += (read_positions.size - keep.sum())
            read_positions = read_positions[keep]
            read_positions[read_positions > (right - adjustment)] -= remove_size

        return read_positions, new_genome_length


    def filter_reads_phase1(self, read_positions, genome_length):
        """A coarse-grained genomewide filter that removes reads in
        ultra-high or ultra-low coverage regions.

        Parameters
        ----------
            read_positions : np.array of int
                The coordinate of the start position of each read
            genome_length : int

        Returns
        -------
            read_positions : np.array of int
                The starting position of each read in an unfiltered region,
                adjusted for the positions removed.
            new_genome_length : float
                The length of the genome with filtered regions removed
        """
        read_positions = np.array(read_positions)
        binned_reads = self.bin_reads_10Kb(read_positions, genome_length)
        median = np.median(binned_reads)
        lower_bound, upper_bound = self.compute_genomewide_bounds(median)

        rolling_counts, endpoints = self.compute_rolling_sum(read_positions, genome_length)

        # find regions to remove
        remove_start = None
        remove_end = None
        remove_regions = []
        for i,c in enumerate(rolling_counts):
            left, right = endpoints[i]

            if (c > upper_bound or c < lower_bound) and remove_start is None:
                remove_start = left

            elif c < upper_bound and c > lower_bound and remove_start is not None:
                remove_end = left # want endpoint of previous bin
                remove_regions.append( (remove_start, remove_end) )
                remove_start = None
                remove_end = None

        # get endpoint
        if remove_start is not None:
            remove_end = genome_length
            remove_regions.append( (remove_start, remove_end) )

        read_positions, new_genome_length = self.remove_reads_by_region(read_positions, genome_length, remove_regions)
        
        binned_reads = self.bin_reads_10Kb(read_positions, genome_length)
        return read_positions, new_genome_length


    def filter_reads_phase2(self, read_positions, genome_length):
        """A fine-grained filter that removes reads in localized regions
        with too-high or too-low coverage. For each 10Kb, looks 6.25% of
        the genome length ahead and 6.25% of the genome length behind.

        Parameters
        ----------
            read_positions : np.array of int
                The coordinate of the start position of each read
            genome_length : int
                The length of the reference genome

        Returns
        -------
            read_positions : np.array of int
                The starting position of each read in an unfiltered region,
                adjusted for the positions removed.
            new_genome_length : float
                The length of the genome with filtered regions removed
        """

        # rolling counts in 10Kb bins, shifting by 1000bp
        rolling_counts, endpoints = self.compute_rolling_sum(read_positions, genome_length)
        nbins = rolling_counts.size

        # how many 10 Kb bins to use
        window = np.max((4,math.ceil(0.0625*genome_length / 10000)))
        window = int(window)

        # find regions to remove
        remove_start = None
        remove_end = None
        remove_regions = []

        bounds = {}

        for i,c in enumerate(rolling_counts):
            # endpoints
            left, right = endpoints[i]

            # take bins in window to the left and right
            left_bins = [(i - 10*j) for j in range(window)]
            right_bins = [(i + 10*j) % nbins for j in range(window)]
            bins = np.concatenate((rolling_counts[left_bins], rolling_counts[right_bins]))
            median = np.median(bins)

            # computing these bounds is expensive, so
            # let's not compute them more than is needed
            if median in bounds:
                lower_bound, upper_bound = bounds[median]
            else:
                lower_bound, upper_bound = self.compute_bin_bounds(median)
                bounds[median] = (lower_bound, upper_bound)


            if (c > upper_bound or c < lower_bound) and remove_start is None:
                remove_start = left

            elif (c < upper_bound and c > lower_bound) and remove_start is not None:
                remove_end = left # want endpoint of previous bin
                remove_regions.append( (remove_start, remove_end) )
                remove_start = None
                remove_end = None

        # get endpoint
        if remove_start is not None:
            remove_end = genome_length
            remove_regions.append((remove_start, remove_end))


        read_positions, new_genome_length = self.remove_reads_by_region(read_positions, genome_length, remove_regions)
        return read_positions, new_genome_length


    def bin_reads_10Kb(self, read_positions, genome_length):
        """Aggregate reads into 10Kb windows and compute the read count in each window.
        
        Parameters
        ----------
            read_positions : np.array of int
                The coordinate of the start position of each read
            genome_length : int
                The length of the reference genome

        Returns
        -------
            bin_counts : np.array
                An array of the read count in each bin
        """
        bin_size = 10000
        nbins = int(math.ceil(genome_length / bin_size))
        bin_counts = np.zeros(nbins)

        for r in read_positions:
            rbin = int(math.floor(nbins * r / genome_length))
            bin_counts[rbin] += 1

        return bin_counts



class CoPTRRefEstimate:
    """Data structure to store CoPTRRef estimates.

    Parameters
    ----------
        bam_file : str
            The bam file from which the estimate was obtained
        genome_id : str
            The reference genome id for the estimate
        sample_id : str
            The sample id for the estimate
        estimate : float or np.nan
            The estimated PTR. This will be set to np.nan if the sample
            did not pass filtering steps
        nreads : int
            The number of reads remaining after filtering
        cov_frac : float
            The fraction of nonzero bins in 10Kb windows.
    """

    def __init__(self, bam_file, genome_id, sample_id, estimate, nreads, cov_frac):
        self.bam_file = bam_file
        self.genome_id = genome_id
        self.sample_id = sample_id
        self.estimate = estimate
        self.nreads = nreads
        self.cov_frac = cov_frac


    def __str__(self):
        return "CoPTRRefEstimate(bam_file={}, genome_id={}, sample_id={}, estimate={:.3f}, nreads={}, cov_frac={}".format(
            self.bam_file, self.genome_id, self.sample_id, self.estimate, self.nreads, self.cov_frac
            )

    def __repr__(self):
        return self.__str__()



class CoPTRRef:
    """CoPTR estimator for complete reference genomes.

    Parameters
    ----------
        min_reads : float
            Minumum read count after filtering to estimate PTR
        frac_nonzero : float
            Fraction of nonzero 10Kb bins required to estimate PTR
    """

    def __init__(self, min_reads, min_cov):
        self.min_reads = min_reads
        self.min_cov = min_cov


    def compute_log_likelihood(self, ori_loc, ter_loc, log2_ptr, read_locs):
        """Model log-likelihood for one sample.

        Parameters
        ------
            ori_loc : float
                Position of replication origin in the interval [0, 1]
            ter_loc : float
                Position of terminus location in the interval [0, 1]
            log2_ptr : float
                The log base 2 peak-to-trough ratio
            read_locs : np.array of float
                A 1D numpy array giving read positions in [0, 1]
        """
        if ori_loc > 1 or ori_loc < 0 or ter_loc > 1 or ter_loc < 0 or log2_ptr > np.log2(4) or log2_ptr <= np.log2(1):
            return -np.inf

        x1 = np.min([ori_loc, ter_loc])
        x2 = np.max([ori_loc, ter_loc])

        dist = np.abs(x2-x1)
        if dist < 0.45 or dist > 0.55:
            return -np.inf

        # compute log2_p(x_i), log2_p(x_t) so that the density normalizes
        alpha = log2_ptr / (ori_loc - ter_loc)
        if ori_loc < ter_loc:
            c1 = np.log2(np.log(2)) - \
                  np.log2((1/alpha)*(2**(alpha*x1) + 2**(alpha*(x2 - x1)) - 2 - 2**(-alpha*(1-x2) - log2_ptr) + 2**(-log2_ptr)))
            c2 = c1 - log2_ptr
        else:
            c1 = np.log2(np.log(2)) - \
                 np.log2((1/alpha)*(2**(alpha*x1) + 2**(alpha*(x2 - x1)) - 2 - 2**(-alpha*(1-x2) + log2_ptr) + 2**(log2_ptr)))
            c2 = c1 + log2_ptr

        log_prob = 0
        if read_locs.size == 1:
            left_part = read_locs[read_locs <= x1]
            middle_part = read_locs[np.logical_and(read_locs > x1, read_locs < x2)]
            right_part = read_locs[read_locs >= x2]
        else:
            left_x = np.searchsorted(read_locs, x1)
            right_x = np.searchsorted(read_locs, x2)
            left_part = read_locs[:left_x]
            middle_part = read_locs[left_x:right_x]
            right_part = read_locs[right_x:]
        log_prob += (-alpha*(left_part - x1)   + c1).sum()
        log_prob += ( alpha*(middle_part - x1) + c1).sum()
        log_prob += (-alpha*(right_part - x2)  + c2).sum()

        return log_prob


    def estimate_ptr(self, read_positions, ref_genome_len):
        """Estimate the PTR for a single sample.

        Parameters
        ------
            read_positions : np.array
                A 1D numpy array giving the starting coordinates of each read
            ref_genome_length : int
                The length of the refernce genome

        Returns
        -------
            log2_ptr : float
                The estimated log base 2 PTR
            ori : float
                The estimated origin in position in [0, 1]
            ter : float
                The estimated replication terminus position in [0, 1]
            log_lk : float
                The model log likelihood
        """
        rf = ReadFilterRef(self.min_reads, self.min_cov)
        read_positions, ref_genome_len, qc_result = rf.filter_reads(read_positions, ref_genome_len)
        if not qc_result.passed_qc:
            return np.nan, np.nan, np.nan, np.nan, qc_result

        read_locs = np.sort(np.array(read_positions) / ref_genome_len)

        # handled by multiple initial conditions
        np.seterr(invalid="ignore")
        f = lambda x: -self.compute_log_likelihood(x[0], x[1], x[2], read_locs)
        # initial ori estimate, initial ptr estimate
        best_x = 0
        best_f = np.inf
        nreads = len(read_locs)

        # try multiple initial conditions
        log_ptrs = [0.2*(i+1) for i in range(5)]
        oris = [0.2*i for i in range(5)]
        np.random.shuffle(log_ptrs)
        np.random.shuffle(oris)
        xs = []
        for log_ptr in log_ptrs:
            for ori in oris:
                ter = (0.5 + ori) % 1
                xs.append([ori, ter, log_ptr])
        xs = np.array(xs)

        for x0 in xs:
            res = scipy.optimize.minimize(f, x0, method="SLSQP")

            if not np.any(np.isnan(res.x)) and res.fun < best_f:
                best_x = res.x
                best_f = res.fun


        np.seterr(invalid="warn")
        if best_f == np.inf:
            print_warning("CoPTRRef", "PTR optimization failed")
            return np.nan, np.nan, np.nan, np.nan
        else:
            ori = np.clip(best_x[0], 0, 1)
            ter = np.clip(best_x[1], 0, 1)
            log2_ptr = np.clip(best_x[2], 0, np.inf)
            return log2_ptr, ori, ter, -best_f, qc_result


    def _parallel_helper(self, x):
        read_positions = x[0]
        genome_length = x[1]
        bam_file = x[2]
        genome_id = x[3]
        sample_id = x[4]
        log2_ptr, ori, ter, log_lk, qc_result = self.estimate_ptr(read_positions, genome_length)
        est = CoPTRRefEstimate(bam_file, genome_id, sample_id, log2_ptr, qc_result.nreads, qc_result.frac_nonzero)
        return est
        


def estimate_ptrs_coptr_ref(coverage_maps, min_reads, min_cov, threads):
    """Estimate Peak-to-Trough ratios across samples.

    Parameters
    ----------
        coverage_maps : dict[str -> list[CoverageMap]]
            A dictionary with key reference genome id and value
            a list of coverage maps for that reference genome
        min_reads : float
            Minumum read count after filtering to estimate PTR
        frac_nonzero : float
            Fraction of nonzero 10Kb bins required to estimate PTR
        threads : int
            Number of threads for parallel computation

    Returns
    -------
        coptr_ref_estimates : dict[str -> list[CoPTRefEstimate]]
            A dictionary with key reference genome id and value
            a list of CoPTRRefEstimate for that reference genome.
    """
    # workers for multiprocessing
    pool = mp.Pool(threads)
    coptr_ref_estimates = {}
    coptr_ref = CoPTRRef(min_reads, min_cov)

    for ref_id in coverage_maps:
        print_info("CoPTRRef", "estimating PTRs for {}".format(ref_id))

        if threads > 1:
            coptr_ref_estimates[ref_id] =  pool.map(coptr_ref._parallel_helper, [ (sample.read_positions, sample.length, sample.bam_file, sample.genome_id, sample.sample_id) for sample in coverage_maps[ref_id]])
        else:
            coptr_ref_estimates[ref_id] = []
            for sample in coverage_maps[ref_id]:
                log2_ptr, ori, ter, log_lk, qc_result = coptr_ref.estimate_ptr(sample.read_positions, sample.length)
                est = CoPTRRefEstimate(sample.bam_file, sample.genome_id, sample.sample_id, log2_ptr, qc_result.nreads, qc_result.frac_nonzero)
                coptr_ref_estimates[ref_id].append(est)
    return coptr_ref_estimates