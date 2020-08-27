"""
coptr_ref.py
======================
Estimate peak-to-trough ratios from assemblies.
"""
import math
import multiprocessing as mp
import numpy as np
import scipy.stats
import sys

from src.print import print_info, print_warning, print_error



class CoPTRContigEstimate:
    """Data structure to store CoPTRContig estimates.

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
            If the sample passing filtering steps (i.e. has a PTR estimate)
            this is the number of filtered reads used for the estimate. Otherwise
            it is the number of reads before filtering.
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
        return "CoPTRContigEstimate(bam_file={}, genome_id={}, sample_id={}, estimate={:.3f}, nreads={}, cov_frac={}".format(
            self.bam_file, self.genome_id, self.sample_id, self.estimate, self.nreads, self.cov_frac
            )

    def __repr__(self):
        return self.__str__()


class CoPTRContig:
    """Estimate PTRs from draft assemblies.
    
    Parameters
    ----------
        min_reads : float
            Minumum read count after filtering to estimate PTR
        min_samples : float
            Minimum number of passing samples required to estimate a PTR
    """

    def __init__(self, min_reads, min_samples):
        self.min_reads = min_reads
        self.min_samples = min_samples


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
        pois_lower = scipy.stats.poisson(0.5*median_read_count)
        pois_upper = scipy.stats.poisson(2*median_read_count)
        for i in range(2, 17):
            m = 0.5*i
            lower_bound = (1./m)*median_read_count
            upper_bound = m*median_read_count
            if pois_upper.cdf(upper_bound) > 0.99 and (1 - pois_lower.cdf(lower_bound)) > 0.99:
                return lower_bound, upper_bound
        return lower_bound, upper_bound


    def construct_coverage_matrix(self, coverage_maps):
        """Aggregate coverage maps into a bin by sample matrix.

        Parameters
        ----------
            coverage_maps : list[CoverageMapContig]
                A list of coverage maps per sample.

        Returns
        -------
            A : 2D np.array of float
                The resulting coverage matrix.
        """
        A = []
        contig_ids = coverage_maps[0].contig_ids
        for cm in coverage_maps:

            binned_reads = []
            for contig_id in sorted(contig_ids):

                if contig_id not in contig_ids:
                    raise print_error("CoPTRContig", "missing contig {} from {}".format(contig_id, cm.bam_file))

                length = cm.get_length(contig_id)
                if length < 11000:
                    continue
                else:
                    binned_reads.append(cm.bin_reads_10Kb(contig_id))
            
            A.append(np.concatenate(binned_reads))

        # bin by sample
        A = np.array(A).T
        return A


    def estimate_slope_intercept(self, bin_probs):
        """Estimate the slope along reordering bins using least squares.

        Parameters
        ----------
            bin_probs : np.array
                The observed probability that a read lands in a 10Kb bin.

        Returns
        -------
            m : float
                The slope
            b : float
                The intercept
        """
        rng = np.linspace(0, 1, bin_probs.size)
        y = bin_probs
        x = np.ones((y.size, 2))
        x[:,0] = rng

        sol, residuals, rank, s = np.linalg.lstsq(x, y, rcond=None)
        m,b = sol

        return m,b


    def estimate_ptrs(self, coverage_maps, return_verbose=False):
        """Estimate PTRs across multiple samples of the same reference genome.
        
        Parameters
        ----------
            coverage_maps : list[CoverageMapContig]
                The coverage maps for the assembly
            min_reads : int
                Minimum number of reads after filtering to estimate a PTR
            min_samples : int
                Minimum number of samples required to reorder bins using PCA.
                No estimate will be produced if the number of passing samples is
                less than min_samples.
            return_verbose : bool
                If true, returns estimated slopes, intercepts, and binned
                read counts in addition to estimates

        Returns
        -------
            estimates : list[CoPTRContigEstimate]
                A list of estimated PTRs
            parameters : list[tuple(2)]
                If return_verbose is true, this is a list of tuples with
                first coordinate the stimated slope, and second coordinate
                the intercept
            binned_counts list[np.array]
                If return_verbose is true, this is a list of the binned 10Kb
                read counts after filtering steps.
        """
        estimates = []
        parameters = []
        binned_counts = []

        passing_coverage_maps = []
        for cm in coverage_maps:
            if cm.passed_qc():
                passing_coverage_maps.append(cm)
            # return nan for estimate
            else:
                bam_file = cm.bam_file
                genome_id = cm.genome_id
                estimate = CoPTRContigEstimate(cm.bam_file, cm.genome_id, cm.sample_id, np.nan, cm.reads, cm.frac_nonzero)
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])

        passing_coverage_maps.sort(key=lambda cm: cm.bam_file)

        # bin by sample
        A = self.construct_coverage_matrix(passing_coverage_maps)

        # remove 10Kb bins without reads
        zero_rows = np.sum(A == 0, axis=1) > 0
        A = A[np.logical_not(zero_rows), :]

        # filter out low coverage and high coverage bins,
        # remove samples with too few reads
        A_filtered = []
        bam_files = []
        genome_ids = []
        sample_ids = []

        tmp = []
        for i,col in enumerate(A.T):
            median = np.median(col)
            lower_bound, upper_bound = self.compute_genomewide_bounds(median)
            col[np.logical_or(col < lower_bound, col > upper_bound)] = 0
            
            if col.sum() > self.min_reads:
                A_filtered.append(col)
                bam_files.append(passing_coverage_maps[i].bam_file)
                genome_ids.append(passing_coverage_maps[i].genome_id)
                sample_ids.append(passing_coverage_maps[i].sample_id)
                tmp.append(passing_coverage_maps[i])
            else:
                bam_file = passing_coverage_maps[i].bam_file
                genome_id = passing_coverage_maps[i].genome_id
                sample_id = passing_coverage_maps[i].sample_id
                reads = passing_coverage_maps[i].reads
                frac_nonzero = passing_coverage_maps[i].frac_nonzero
                m = np.nan
                estimate = CoPTRContigEstimate(bam_file, genome_id, sample_id, np.nan, col.sum(), frac_nonzero)
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])

        if len(A_filtered) < self.min_samples:
            for cm in tmp:
                bam_file = cm.bam_file
                genome_id = cm.genome_id
                estimate = CoPTRContigEstimate(cm.bam_file, cm.genome_id, cm.sample_id, np.nan, cm.reads, cm.frac_nonzero)
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])
            return estimates
        else:
            tmp = None

        A = np.array(A_filtered).T

        # bins to be filtered are set to 0, so let's remove zeros again
        zero_rows = np.sum(A == 0, axis=1) > 0
        A = A[np.logical_not(zero_rows), :]

        # normalize so the bin counts are now probabilities
        A /= A.sum(axis=0,keepdims=True)
        A = np.log2(A)

        # project bins onto first PC
        u,s,vh = np.linalg.svd(A)
        pc1 = A.dot(vh[0])

        # reorder bins by PC
        sorted_idx = np.argsort(pc1)
        A = A[sorted_idx,:]

        for j,col in enumerate(A.T):
            lb = int(0.05*col.size)
            ub = int(0.95*col.size)
            m,b = self.estimate_slope_intercept(col[lb:ub])
            parameters.append((m, b))

            # the sign of the slope depends on the orientation,
            # so let's look at positive slopes only
            m = np.abs(m)
            estimate = CoPTRContigEstimate(bam_files[j], genome_ids[j], sample_ids[j], m, col.sum(), (col == 0).sum() / col.size)
            estimates.append(estimate)

            binned_counts.append(np.power(2,col))

        if return_verbose:
            return estimates, parameters, binned_counts
        else:
            return estimates

    def _parallel_helper(self, x):
        ref_genome = x[0]
        coverage_maps = x[1]

        return (ref_genome, self.estimate_ptrs(coverage_maps))



def estimate_ptrs_coptr_contig(coverage_maps, min_reads, min_samples, threads):
    """Estimate Peak-to-Trough ratios across samples.

    Parameters
    ----------
        coverage_maps : dict[str -> list[CoverageMap]]
            A dictionary with key reference genome id and value
            a list of coverage maps for that reference genome
        min_reads : float
            Minumum read count after filtering to estimate PTR
        min_samples : float
            Minimum number of passing samples required to estimate a PTR
        threads : int
            Number of threads for parallel computation

    Returns
    -------
        coptr_contig_estimates : dict[str -> list[CoPTContigEstimate]]
            A dictionary with key reference genome id and value
            a list of CoPTRContigEstimate for that reference genome.
    """
    coptr_contig_estimates = {}
    coptr_contig = CoPTRContig(min_reads, min_samples)
    print_info("CoPTRContig", "estimating PTRs across genomes")

    if threads == 1:
        for ref_id in coverage_maps:
            coptr_contig_estimates[ref_id] = coptr_contig.estimate_ptrs(coverage_maps[ref_id])
    else:
        # workers for multiprocessing
        pool = mp.Pool(threads)
        flat_coverage_maps = [(ref_id, coverage_maps[ref_id]) for ref_id in coverage_maps]
        flat_results = pool.map(coptr_contig._parallel_helper, flat_coverage_maps)
        coptr_contig_estimates = dict(flat_results)

    return coptr_contig_estimates