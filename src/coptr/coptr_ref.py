"""
coptr_ref.py
======================
Estimate peak-to-trough ratios using complete reference genomes.
"""

"""
This file is part of CoPTR.

CoPTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CoPTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoPTR.  If not, see <https://www.gnu.org/licenses/>.
"""

import logging
import math
import multiprocessing as mp
import os.path
import pickle as pkl
import struct
import sys

import numpy as np
import scipy.optimize
import scipy.stats


logger = logging.getLogger(__name__)


class QCResult:
    """Data structure to store results of quality checking.

    Parameters
    ----------
        frac_nonzero : float
            Fraction of nonzero bins in 10Kb windows
        nreads : int
            Number of reads after filtering
        frac_removed : float
            Proportion of genome removed during filtering
        passed_qc : bool
            Flag indicating if the sample passed quality thresholds
    """

    def __init__(self, frac_nonzero, nreads, frac_removed, passed_qc):
        self.frac_nonzero = frac_nonzero
        self.nreads = nreads
        self.frac_removed = frac_removed
        self.passed_qc = passed_qc

    def __str__(self):
        return "QCResult(frac_nonzero={}, nreads={}, frac_removed={}, passed_qc={})".format(
            self.frac_nonzero, self.nreads, self.frac_removed, self.passed_qc
        )

    def __repr__(self):
        return self.__str__()


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
        # don't filter if there are too few reads
        # just return that the sample failed QC
        if len(read_positions) < self.min_reads:
            return (
                np.array([]),
                genome_length,
                QCResult(-1, len(read_positions), -1, False),
            )

        bin_size = self.compute_bin_size(genome_length)
        read_positions = np.copy(read_positions)
        filtered_read_positions, filtered_genome_length = self.filter_reads_phase1(
            read_positions, genome_length, bin_size
        )
        # more than 25% of the genome is removed
        frac_removed = 1 - filtered_genome_length / genome_length
        if frac_removed > 0.25:
            return (
                np.array(filtered_read_positions),
                filtered_genome_length,
                QCResult(-1, filtered_read_positions.size, frac_removed, False),
            )

        filtered_read_positions, filtered_genome_length = self.filter_reads_phase2(
            filtered_read_positions, filtered_genome_length, bin_size
        )
        qc_result = self.quality_check(
            filtered_read_positions,
            filtered_genome_length,
            genome_length,
            bin_size,
            self.min_reads,
            self.min_cov,
        )
        return filtered_read_positions, filtered_genome_length, qc_result

    def compute_bin_size(self, genome_length):
        """Compute bin size for read counts.

        Parameters
        ----------
            genome_length : int
                Length of the reference genome

        Returns
        -------
            bin_size : int
                Bin size for read countns
        """
        # we want approximately 500 bins
        target_bins = 500
        bin_size = genome_length / target_bins

        # want a number divisible by 100 for downstream steps
        bin_size = bin_size - (bin_size % 100)

        if bin_size < 1:
            logger.error(
                "Found complete reference genome with <500bp. "
                "This is probably due to a mislabeled contig. "
                "Please check your .genomes file."
            )

        return bin_size

    def quality_check(
        self,
        filtered_read_positions,
        filtered_genome_length,
        original_genome_length,
        bin_size,
        min_reads,
        frac_nonzero,
    ):
        """A basic quality check for required coverage of a genome in a sample.

        Parameters
        ----------
            filtered_read_positions : np.array or list of int
                The starting coordinate of each read along the genome
                after filtering
            filtered_genome_length : float
                The length of the reference genome after filtering
            original_genome_length : float
                The length of the reference genome before filtering
            min_reads : float
                Minumum read count after filtering
            frac_nonzero : float
                Fraction of nonzero 10Kb bins required

        Returns
        -------
            qc_result : QCResult
                Object that stores the result of the quality check
        """
        binned_counts = self.bin_reads(
            filtered_read_positions, filtered_genome_length, bin_size
        )
        nonzero_bins = (binned_counts > 0).sum() / np.max((1, binned_counts.size))
        nreads = binned_counts.sum()

        passed_qc = (
            True if nonzero_bins >= frac_nonzero and nreads >= min_reads else False
        )

        frac_removed = 1 - filtered_genome_length / original_genome_length
        if frac_removed > 0.25:
            passed_qc = False

        return QCResult(nonzero_bins, nreads, frac_removed, passed_qc)

    def compute_genomewide_bounds(self, binned_counts):
        """Compute bounds on read counts in bins.

        Parameters
        ----------
            binned_counts : np.array
                Binned read counts

        Returns
        -------
            lower_bound : float
                A lower bound on read counts in bins
            upper_bound : float
                An upper bound on read counts in bins

        """
        median = np.log2(np.median(binned_counts[binned_counts != 0]))

        std = np.std(np.log2(binned_counts[binned_counts != 0]))
        std = np.max((1, std))

        norm = scipy.stats.norm()
        m = norm.ppf(0.95)

        lower_bound = np.power(2, median - m * std)
        upper_bound = np.power(2, median + m * std)

        return lower_bound, upper_bound

    def compute_bin_bounds(self, binned_counts):
        """Compute bounds on read counts in a smaller window along the genome.

        Parameters
        ----------
            binned_counts : np.array
                Binned read counts

        Returns
        -------
            lower_bound : float
                A lower bound on read counts in bins
            upper_bound : float
                An upper bound on read counts in bins
        """
        median = np.log2(np.median(binned_counts[binned_counts != 0]))
        std = np.std(np.log2(binned_counts[binned_counts != 0]))
        std = np.max((1, std))

        norm = scipy.stats.norm()
        m = norm.ppf(0.975)

        lower_bound = np.power(2, median - m * std)
        upper_bound = np.power(2, median + m * std)

        return lower_bound, upper_bound

    def compute_rolling_sum(self, read_positions, genome_length, bin_size):
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
        step = math.ceil(bin_size / 100)
        s = read_positions.size

        sorted_reads = np.sort(read_positions)
        endpoints = np.array([(i, i + bin_size) for i in range(0, genome_length, step)])

        left_endpoints = [e[0] for e in endpoints]
        right_endpoints = [e[1] for e in endpoints]
        left_idx = np.searchsorted(sorted_reads, left_endpoints, side="right")
        right_idx = np.searchsorted(sorted_reads, right_endpoints, side="right")
        rolling_counts = right_idx - left_idx

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
        for left, right in regions:
            remove_size = int(right - left)
            new_genome_length -= remove_size

            keep = np.logical_or(
                read_positions < (left - adjustment),
                read_positions >= (right - adjustment),
            )
            adjustment += remove_size
            filter_count += read_positions.size - keep.sum()
            read_positions = read_positions[keep]
            read_positions[read_positions > (right - adjustment)] -= remove_size

        return read_positions, new_genome_length

    def filter_reads_phase1(self, read_positions, genome_length, bin_size):
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
        binned_reads = self.bin_reads(read_positions, genome_length, bin_size)
        lower_bound, upper_bound = self.compute_genomewide_bounds(binned_reads)

        rolling_counts, endpoints = self.compute_rolling_sum(
            read_positions, genome_length, bin_size
        )

        # find regions to remove
        remove_start = None
        remove_end = None
        remove_regions = []
        for i, c in enumerate(rolling_counts):
            left, right = endpoints[i]

            if (c > upper_bound or c < lower_bound) and remove_start is None:
                remove_start = left

            elif c < upper_bound and c > lower_bound and remove_start is not None:
                remove_end = left  # want endpoint of previous bin
                remove_regions.append((remove_start, remove_end))
                remove_start = None
                remove_end = None

        # get endpoint
        if remove_start is not None:
            remove_end = genome_length
            remove_regions.append((remove_start, remove_end))

        read_positions, new_genome_length = self.remove_reads_by_region(
            read_positions, genome_length, remove_regions
        )

        binned_reads = self.bin_reads(read_positions, genome_length, bin_size)
        return read_positions, new_genome_length

    def filter_reads_phase2(self, read_positions, genome_length, bin_size):
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
        rolling_counts, endpoints = self.compute_rolling_sum(
            read_positions, genome_length, bin_size
        )
        nbins = rolling_counts.size

        # how many bins to use
        window = np.max((4, math.ceil(0.0625 * genome_length / bin_size)))
        window = int(window)

        # find regions to remove
        remove_start = None
        remove_end = None
        remove_regions = []

        bounds = {}

        for i, c in enumerate(rolling_counts):
            # endpoints
            left, right = endpoints[i]

            # take bins in window to the left and right
            left_bins = [(i - 100 * j) for j in range(window)]
            right_bins = [(i + 100 * j) % nbins for j in range(window)]
            bins = np.concatenate(
                (rolling_counts[left_bins], rolling_counts[right_bins])
            )
            median = np.median(bins)

            # computing these bounds is expensive, so
            # let's not compute them more than is needed
            if median in bounds:
                lower_bound, upper_bound = bounds[median]
            else:
                lower_bound, upper_bound = self.compute_bin_bounds(bins)
                bounds[median] = (lower_bound, upper_bound)

            if (c > upper_bound or c < lower_bound) and remove_start is None:
                remove_start = left

            elif (c < upper_bound and c > lower_bound) and remove_start is not None:
                remove_end = left  # want endpoint of previous bin
                remove_regions.append((remove_start, remove_end))
                remove_start = None
                remove_end = None

        # get endpoint
        if remove_start is not None:
            remove_end = genome_length
            remove_regions.append((remove_start, remove_end))

        read_positions, new_genome_length = self.remove_reads_by_region(
            read_positions, genome_length, remove_regions
        )
        return read_positions, new_genome_length

    def bin_reads(self, read_positions, genome_length, bin_size=10000):
        """Aggregate reads into bin_size windows and compute the read count in each window.

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
            The estimated log2(PTR). This will be set to np.nan if the sample
            did not pass filtering steps
        ori_estimate : float
            Estimated replication origin position in interval [0, 1]
        ter_estimate : float
            Estimated replication terminus position in interval [0, 1]
        nreads : int
            The number of reads remaining after filtering
        cov_frac : float
            The fraction of nonzero bins in 10Kb windows.
    """

    def __init__(
        self,
        bam_file,
        genome_id,
        sample_id,
        estimate,
        ori_estimate,
        ter_estimate,
        nreads,
        cov_frac,
    ):
        self.bam_file = bam_file
        self.genome_id = genome_id
        self.sample_id = sample_id
        self.estimate = estimate
        self.ori_estimate = ori_estimate
        self.ter_estimate = ter_estimate
        self.nreads = nreads
        self.cov_frac = cov_frac

    def __str__(self):
        return "CoPTRRefEstimate(bam_file={}, genome_id={}, sample_id={}, estimate={:.3f}, nreads={}, cov_frac={})".format(
            self.bam_file,
            self.genome_id,
            self.sample_id,
            self.estimate,
            self.nreads,
            self.cov_frac,
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
        if (
            ori_loc > 1
            or ori_loc < 0
            or ter_loc > 1
            or ter_loc < 0
            or log2_ptr > np.log2(16)
            or log2_ptr <= np.log2(1)
        ):
            return -np.inf

        assert read_locs.size == 1 or np.all(
            read_locs[:-1] <= read_locs[1:]
        ), "reads must be sorted"

        x1 = np.min([ori_loc, ter_loc])
        x2 = np.max([ori_loc, ter_loc])

        dist = np.abs(x2 - x1)
        if dist < 0.45 or dist > 0.55:
            return -np.inf

        # compute log2_p(x_i), log2_p(x_t) so that the density normalizes
        alpha = log2_ptr / (ori_loc - ter_loc)
        if ori_loc < ter_loc:
            c1 = np.log2(np.log(2)) - np.log2(
                (1 / alpha)
                * (
                    2 ** (alpha * x1)
                    + 2 ** (alpha * (x2 - x1))
                    - 2
                    - 2 ** (-alpha * (1 - x2) - log2_ptr)
                    + 2 ** (-log2_ptr)
                )
            )
            c2 = c1 - log2_ptr
        else:
            c1 = np.log2(np.log(2)) - np.log2(
                (1 / alpha)
                * (
                    2 ** (alpha * x1)
                    + 2 ** (alpha * (x2 - x1))
                    - 2
                    - 2 ** (-alpha * (1 - x2) + log2_ptr)
                    + 2 ** (log2_ptr)
                )
            )
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
        log_prob += (-alpha * (left_part - x1) + c1).sum()
        log_prob += (alpha * (middle_part - x1) + c1).sum()
        log_prob += (-alpha * (right_part - x2) + c2).sum()

        return log_prob

    def compute_multisample_log_likelihood(
        self, ori_loc, ter_loc, log2_ptrs, read_locs_list
    ):
        log_lk = 0
        for log2_ptr, read_locs in zip(log2_ptrs, read_locs_list):
            log_lk += self.compute_log_likelihood(ori_loc, ter_loc, log2_ptr, read_locs)
        return log_lk

    def estimate_ptr(
        self, read_positions, ref_genome_len, filter_reads=True, estimate_terminus=False
    ):
        """Estimate the PTR for a single sample.

        Parameters
        ------
            read_positions : np.array
                A 1D numpy array giving the starting coordinates of each read
            ref_genome_length : int
                The length of the refernce genome
            filter_reads : bool
                If true, reads a filtered before estimates
            estimate_terminus :
                If true, the replication terminus is estimated in addition to the
                replication origin

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
        if filter_reads:
            rf = ReadFilterRef(self.min_reads, self.min_cov)
            read_positions, ref_genome_len, qc_result = rf.filter_reads(
                read_positions, ref_genome_len
            )
            if not qc_result.passed_qc:
                return np.nan, np.nan, np.nan, np.nan, qc_result

        read_locs = np.sort(np.array(read_positions) / ref_genome_len)

        # handled by multiple initial conditions
        np.seterr(invalid="ignore")

        if estimate_terminus:
            f = lambda x: -self.compute_log_likelihood(x[0], x[1], x[2], read_locs)

            log_ptrs = [0.2 * (i + 1) for i in range(5)]
            oris = [0.2 * i for i in range(5)]
            np.random.shuffle(log_ptrs)
            np.random.shuffle(oris)
            xs = []
            for log_ptr in log_ptrs:
                for ori in oris:
                    ter = (0.5 + ori) % 1
                    xs.append([ori, ter, log_ptr])
            xs = np.array(xs)

        else:
            f = lambda x: -self.compute_log_likelihood(
                x[0], (x[0] + 0.5) % 1, x[1], read_locs
            )

            log_ptrs = [0.2 * (i + 1) for i in range(5)]
            oris = [0.2 * i for i in range(5)]
            np.random.shuffle(log_ptrs)
            np.random.shuffle(oris)
            xs = []
            for log_ptr in log_ptrs:
                for ori in oris:
                    xs.append([ori, log_ptr])
            xs = np.array(xs)

        # initial ori estimate, initial ptr estimate
        best_x = 0
        best_f = np.inf
        nreads = len(read_locs)

        for x0 in xs:
            res = scipy.optimize.minimize(f, x0, method="SLSQP")

            if not np.any(np.isnan(res.x)) and res.fun < best_f:
                best_x = res.x
                best_f = res.fun

        if estimate_terminus:
            ori = np.clip(best_x[0], 0, 1)
            ter = np.clip(best_x[1], 0, 1)
            log2_ptr = np.clip(best_x[2], 0, np.inf)
        else:
            ori = np.clip(best_x[0], 0, 1)
            ter = (ori + 0.5) % 1
            log2_ptr = np.clip(best_x[1], 0, np.inf)

        np.seterr(invalid="warn")
        if best_f == np.inf:
            logger.error("PTR optimization failed.")
            return np.nan, np.nan, np.nan, np.nan
        elif not filter_reads:
            return log2_ptr, ori, ter, -best_f
        else:
            return log2_ptr, ori, ter, -best_f, qc_result

    def update_ori_ter(self, log2_ptrs, read_locs_list):
        """Compute maximum likelihood ori estimates across samples
        given estimated log2(PTR)s across samples.

        Parameters
        ----------
            log2_ptrs : list[float]
                List of estimated log2(PTR)s across samples
            read_locs_list list[np.array]
                List of read positions in [0, 1] across samples

        Returns
        -------
            ori : float
                The estimated replication origin in [0, 1]
            ter : float
                The estimated replication terminus in [0, 1]

        """
        np.seterr(invalid="ignore")
        xs = [0.1 * i for i in range(10)]
        f = lambda x: -self.compute_multisample_log_likelihood(
            x[0], (x[0] + 0.5) % 1, log2_ptrs, read_locs_list
        )

        best_x = None
        best_f = np.inf
        for x0 in xs:
            res = scipy.optimize.minimize(f, x0, method="SLSQP")

            if not np.any(np.isnan(res.x)) and res.fun < best_f:
                best_x = res.x
                best_f = res.fun

        np.seterr(invalid="warn")
        return best_x, (best_x + 0.5) % 1

    def update_ptrs(self, ori, ter, prv_log2_ptrs, read_locs_list):
        """Compute maximum likelihood PTR estimates across samples
        given an estimated replication origin and terminus.

        Parameters
        ----------
            ori : float
                The estimated origin in position in [0, 1]
            ter : float
                The estimated replication terminus position in [0, 1]
            prv_log2_ptrs : list[float]
                Estimated log2(PTR) across samples
            read_locs_list list[np.array]
                Read positions in [0, 1] across samples

        Returns
        -------
            log2_ptrs : list[float]
                The updated log2(PTR) estimates
        """

        np.seterr(invalid="ignore")
        log2_ptrs = []
        for i, read_locs in enumerate(read_locs_list):
            f = lambda x: -self.compute_log_likelihood(ori, ter, x, read_locs)
            res = scipy.optimize.minimize(f, prv_log2_ptrs[i], method="SLSQP")
            log2_ptrs.append(res.x[0])
        np.seterr(invalid="warn")
        return log2_ptrs

    def estimate_ptrs(self, coverage_maps):
        """Compute maximum likelihood PTR estimates across samples.

        Parameters
        ----------
            coverage_maps : list[CoverageMapRef]
                A list of coverage maps per sample

        Returns
        -------
            estimates : list[CoPTRefEstimate]
                A list of estimates per sample

        """
        rf = ReadFilterRef(self.min_reads, self.min_cov)
        estimates = []
        sample_ids = []
        read_positions_list = []
        lengths = []
        genome_id = coverage_maps[0].genome_id
        for cm in coverage_maps:
            read_positions, ref_genome_len, qc_result = rf.filter_reads(
                cm.read_positions, cm.length
            )
            if qc_result.passed_qc:
                read_positions_list.append(np.array(read_positions))
                lengths.append(ref_genome_len)
                sample_ids.append(cm.sample_id)

            estimates.append(
                CoPTRRefEstimate(
                    cm.bam_file,
                    cm.genome_id,
                    cm.sample_id,
                    np.nan,
                    np.nan,
                    np.nan,
                    qc_result.nreads,
                    qc_result.frac_nonzero,
                )
            )

        if len(sample_ids) == 0:
            return estimates

        logger.info("Running %s.", genome_id)

        # first, compute individiual log2_ptr, ori, ter estimates
        log2_ptrs = []
        read_locs_list = []
        for read_positions, length in zip(read_positions_list, lengths):
            log2_ptr, ori, ter, f = self.estimate_ptr(
                read_positions, length, filter_reads=False
            )
            log2_ptrs.append(log2_ptr)
            read_locs_list.append(np.sort(read_positions / length))

        if len(read_locs_list) > 1:
            # update replication origin
            ori, ter = self.update_ori_ter(log2_ptrs, read_locs_list)
            # update ptrs
            log2_ptrs = self.update_ptrs(ori, ter, log2_ptrs, read_locs_list)

        n = 0
        for i, cm in enumerate(coverage_maps):
            if cm.sample_id in sample_ids:
                # sanity check that the ptr estimate is good
                if np.abs(log2_ptrs[n] - 4) < 1e-4:
                    estimates[i].ori_estimate = ori
                    estimates[i].ter_estimate = ter
                    estimates[i].estimate = np.nan
                else:
                    estimates[i].ori_estimate = ori
                    estimates[i].ter_estimate = ter
                    estimates[i].estimate = log2_ptrs[n]
                n += 1

        logger.info("Finished %s.", genome_id)

        return estimates

    def _parallel_helper(self, x):
        ref_genome = x[0]
        coverage_maps = x[1]

        return (ref_genome, self.estimate_ptrs(coverage_maps))


def plot_fit(
    coptr_ref_est, read_positions, genome_length, min_reads, min_cov, plot_folder
):
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        logger.critical("Please install matplotlib to enable plotting.")
        raise

    coptr = CoPTRRef(min_reads, min_cov)
    rf = ReadFilterRef(min_reads, min_cov)

    rp = np.array(read_positions)
    length = genome_length

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))

    # first plot unfiltered bins
    bin_size = rf.compute_bin_size(genome_length)
    binned_counts = rf.bin_reads(rp, length, bin_size)
    x = np.linspace(0, 1, binned_counts.size)
    median = np.median(binned_counts[binned_counts != 0])
    lower_bound, upper_bound = rf.compute_genomewide_bounds(binned_counts)

    min_y = np.min(binned_counts[binned_counts != 0])
    min_y -= 0.5 * min_y
    max_y = np.max(binned_counts[binned_counts != 0]) + 0.008 * binned_counts.sum()

    # plot bins
    ax[0].scatter(x, binned_counts)
    ax[0].set_yscale("log", base=2)
    ax[0].set_ylim((min_y, max_y))
    y_lim = ax[0].get_ylim()

    # plot upper and lower bounds
    ax[0].plot(x, median * np.ones(binned_counts.size), c="red", linewidth=3)
    ax[0].plot(
        x, lower_bound * np.ones(binned_counts.size), c="red", ls=":", linewidth=3
    )
    ax[0].plot(
        x, upper_bound * np.ones(binned_counts.size), c="red", ls=":", linewidth=3
    )

    ax[0].set_ylabel("Read Count")
    ax[0].set_xlabel("Position Along Genome")

    # now plot filtered bins and model fit
    filtered_reads, filtered_length, qc_result = rf.filter_reads(rp, genome_length)
    binned_counts = rf.bin_reads(filtered_reads, filtered_length, bin_size)

    x = np.linspace(0, 1, binned_counts.size)
    probs = binned_counts / binned_counts.sum()
    ax[1].scatter(x, probs, c="C1")
    ax[1].set_yscale("log", base=2)
    min_y = np.min(probs[probs != 0])
    min_y -= 0.5 * min_y
    max_y = np.max(probs) + 0.008
    ax[1].set_ylim((min_y, max_y))

    # plot fit
    y = [
        coptr.compute_log_likelihood(
            coptr_ref_est.ori_estimate,
            coptr_ref_est.ter_estimate,
            coptr_ref_est.estimate,
            xi,
        )
        for xi in x
    ]
    y = np.power(2, y) * (x[1] - x[0])
    ax[1].plot(x, y, c="black", linewidth=3)
    ax[1].set_ylabel("Density")
    ax[1].set_xlabel("Position Along Genome")
    ax[1].set_title(
        "\nlog2(PTR)={:.3f} (Reads = {})".format(
            coptr_ref_est.estimate, int(binned_counts.sum())
        )
    )

    plt.tight_layout()
    plt.savefig(
        os.path.join(
            plot_folder,
            coptr_ref_est.sample_id + "-" + coptr_ref_est.genome_id + ".pdf",
        )
    )
    plt.close()


def load_coverage_maps(file_handle, ref_id):
    coverage_maps = []
    try:
        while True:
            coverage_maps.append(pkl.load(file_handle))

    except EOFError:
        pass
    return coverage_maps


def estimate_ptrs_coptr_ref(
    ref_genome_ids,
    grouped_coverage_map_folder,
    min_reads,
    min_cov,
    plot_folder=None,
    mem_limit=None,
):
    """Estimate Peak-to-Trough ratios across samples.

    Parameters
    ----------
        ref_genome_ids : set[str]
            A list of genome ids to estimate PTRs.
        grouped_coverage_map_folder : str
            Path to folder containing coverage maps grouped by genome.
        min_reads : float
            Minumum read count after filtering to estimate PTR
        frac_nonzero : float
            Fraction of nonzero 10Kb bins required to estimate PTR
        threads : int
            Number of threads for parallel computation
        plot_folder : str
            If not None, plots of fitted model are generated and saved here
        mem_limit : int
            Limit amount of coverage maps loaded into memory at once.
            Breaks samples per genome into batches. Limit is specified in GB.
            Works only in single-threaded mode.

    Returns
    -------
        coptr_ref_estimates : dict[str -> list[CoPTRefEstimate]]
            A dictionary with key reference genome id and value
            a list of CoPTRRefEstimate for that reference genome.
    """
    logger.info("Checking reference genomes.")
    coptr_ref_estimates = {}
    coptr_ref = CoPTRRef(min_reads, min_cov)

    for ref_id in sorted(ref_genome_ids):
        with open(
            os.path.join(grouped_coverage_map_folder, ref_id + ".cm.pkl"), "rb"
        ) as f:
            coverage_maps = load_coverage_maps(f, ref_id)

        coptr_ref_estimates[ref_id] = coptr_ref.estimate_ptrs(coverage_maps)

        for sample, est in zip(coverage_maps, coptr_ref_estimates[ref_id]):
            # only plot samples with a PTR estimate
            if plot_folder is not None and not np.isnan(est.estimate):
                plot_fit(
                    est,
                    sample.read_positions,
                    sample.length,
                    min_reads,
                    min_cov,
                    plot_folder,
                )

    return coptr_ref_estimates
