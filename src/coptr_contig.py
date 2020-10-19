"""
coptr_contig.py
======================
Estimate peak-to-trough ratios from assemblies.
"""
import math
import multiprocessing as mp
import numpy as np
import os.path
import pickle as pkl
import scipy.stats
import sys

from src.poisson_pca import PoissonPCA
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
        m = norm.ppf(0.925)

        lower_bound = np.power(2, median - m*std)
        upper_bound = np.power(2, median + m*std)

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

        # find all common contig ids
        contig_id_list = [cm.contig_ids for cm in coverage_maps]
        contig_ids = coverage_maps[0].contig_ids
        ref_genome = coverage_maps[0].genome_id

        for cm in coverage_maps:

            binned_reads = []
            for contig_id in sorted(contig_ids):

                if contig_id not in cm.contig_ids:
                    print_error("CoPTRContig", "missing contig {} from {} in {}".format(contig_id, ref_genome, cm.sample_id), quit=True)

                length = cm.get_length(contig_id)
                if length < 11000:
                    continue
                else:
                    binned_reads.append(cm.bin_reads(contig_id))

            row = np.concatenate(binned_reads)
            A.append(row)

        # bin by sample
        A = np.array(A).T
        return A


    def estimate_slope_intercept(self, bin_log_probs):
        """Estimate the slope along reordering bins using least squares.

        Parameters
        ----------
            bin_log_probs : np.array
                The observed probability that a read lands in a bin.

        Returns
        -------
            m : float
                The slope
            b : float
                The intercept
        """
        rng = np.linspace(0, 1, bin_log_probs.size)
        y = bin_log_probs
        x = np.ones((y.size, 2))
        x[:,0] = rng

        sol, residuals, rank, s = np.linalg.lstsq(x, y, rcond=None)
        m,b = sol

        pred = m*rng + b
        residuals = y - pred

        return m,b,residuals


    def compute_log_likelihood(self, log2_ptr, binned_counts):
        """Compute the log likelihood of reordered bin probabilities given a log2(PTR).

        Parameters
        ----------
            log2_ptr : float
                The Log2(PTR)
            binned_counts : np.array
                Binned read counts along the genome

        Returns
        -------
            log_lk : float
                The log likelihood
        """
        # origin and terminus are at positions 0 and 1
        ori_loc = 0
        ter_loc = 1

        x = np.linspace(0, 1, binned_counts.size)
        # compute normalized probabilities
        log_probs = np.array([log2_ptr*xi for xi in x])
        log_norm = np.log2(np.power(2, log_probs).sum())
        log_probs -= log_norm

        log_lk = (binned_counts*log_probs).sum()

        return log_lk


    def estimate_ptr_maximum_likelihood(self, binned_counts):
        """Compute maximum likelihood estimates of log2(PTR) gived observed reads.

        Parameters
        ----------
            binned_counts : np.array
                Binned read counts for a genome

        Returns
        -------
            log2_ptr : float
                The maximum likelihood log2(PTR) estimate
            intercept : float
                Probability a reads is in the leftmost bin
            log_lk : float
                The log likelihood

        """
        f = lambda x: -self.compute_log_likelihood(x, binned_counts)
        x0 = 0.5

        res = scipy.optimize.minimize_scalar(f, method="brent")
        log2_ptr = res.x

        x = np.linspace(0, 1, binned_counts.size)
        # compute normalized probabilities
        log_probs = np.array([log2_ptr*xi for xi in x])
        log_norm = np.log2(np.power(log_probs, 2).sum())
        log_probs -= log_norm
        probs = np.power(2, log_probs)
        intercept = probs[0]

        log_lk = -f(log2_ptr)

        # the orientation of the bins is arbitrary,
        # so let's orient with a positive slope
        if log2_ptr < 0:
            intercept = log2_ptr + intercept
            log2_ptr = np.abs(log2_ptr)

        return log2_ptr, intercept, log_lk


    def estimate_ptrs(self, coverage_maps, return_bins=False):
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
            return_bins : bool
                If true, returns estimated slopes, intercepts, and binned
                read counts in addition to estimates

        Returns
        -------
            estimates : list[CoPTRContigEstimate]
                A list of estimated PTRs
            parameters : list[tuple(2)]
                If return_bins is true, this is a list of tuples with
                first coordinate the stimated slope, and second coordinate
                the intercept
            binned_counts list[np.array]
                If return_bins is true, this is a list of the binned 10Kb
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

        if len(passing_coverage_maps) == 0:
            if return_bins:
                return estimates, parameters, binned_counts
            else:
                return estimates

        genome_id = coverage_maps[0].genome_id

        # bin by sample
        A = self.construct_coverage_matrix(passing_coverage_maps)
        nbins = A.shape[0]


        # filter out low coverage and high coverage bins,
        # remove samples with too few reads
        A_filtered = []
        bam_files = []
        genome_ids = []
        sample_ids = []

        tmp = []
        np.set_printoptions(suppress=True)
        for i,col in enumerate(A.T):
            lower_bound, upper_bound = self.compute_genomewide_bounds(col)
            col[np.logical_or(col < lower_bound, col > upper_bound)] = -1

            frac_filtered = (col == -1).sum() / col.size
            
            if col.sum() > self.min_reads and frac_filtered < 0.25:
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
            if return_bins:
                return estimates, parameters, binned_counts
            else:
                return estimates
        else:
            tmp = None

        print_info("CoPTRContig", "running {}".format(genome_id))

        A = np.array(A_filtered).T

        # bins to be filtered are set to -1, so let's remove them
        filter_rows = np.sum(A == -1, axis=1) > 0
        A = A[np.logical_not(filter_rows), :]

        poisson_pca = PoissonPCA()
        W,V = poisson_pca.pca(A,k=1)
        sorted_idx = np.argsort(W.flatten())
        A = A[sorted_idx,:]
        reads = A.sum(axis=0)

        for j,col in enumerate(A.T):
            lb = int(0.05*col.size)
            ub = int(0.95*col.size)
            col = col[lb:ub]

            m,b,log_lk = self.estimate_ptr_maximum_likelihood(col)

            if m < 0:
                flip = True
                b = m + b
                m = np.abs(m)
            else:
                flip = False

            # the sign of the slope depends on the orientation,
            # so let's look at positive slopes only
            m = np.abs(m)
            estimate = CoPTRContigEstimate(bam_files[j], genome_ids[j], sample_ids[j], m, reads[j], col.sum())
            estimates.append(estimate)

            parameters.append((m, b))

            #bc = np.power(2,col)
            bc = col
            if flip:
                bc = np.flip(bc)
            binned_counts.append(bc)


        print_info("CoPTRContig", "finished {}".format(genome_id))

        if return_bins:
            return estimates, parameters, binned_counts
        else:
            return estimates

    def _parallel_helper(self, x):
        ref_genome = x[0]
        coverage_maps = x[1]
        plot_folder = x[2]

        if plot_folder is not None:
            estimates, parameters, reordered_bins = self.estimate_ptrs(coverage_maps, return_bins=True)
            plot_fit(estimates, parameters, coverage_maps, reordered_bins, plot_folder)
        else:
            estimates = self.estimate_ptrs(coverage_maps)

        return (ref_genome, estimates)



def plot_fit(estimates, parameters, coverage_maps, reordered_bins, plot_folder):
    import matplotlib.pyplot as plt

    coptr_contig = CoPTRContig(min_reads=5000, min_samples=5)

    # make sure esimates and coverage maps are in the same order
    order = dict([ (est.sample_id, i) for i,est in enumerate(estimates)])
    tmp = [[] for cm in coverage_maps]
    for cm in coverage_maps:
        idx = order[cm.sample_id]
        tmp[idx] = cm
    coverage_maps = tmp


    nplots = 0
    for estimate,p,cm,rbins in zip(estimates, parameters, coverage_maps, reordered_bins):
        if np.isnan(estimate.estimate):
            continue

        contig_ids = cm.contig_ids
        lengths = []
        binned_counts = []
        for cid in contig_ids:
            length = cm.get_length(cid)
            if length < 11000: continue

            binned_counts += cm.bin_reads(cid).tolist()

        binned_counts = np.array(binned_counts)
        lower_bound, upper_bound = coptr_contig.compute_genomewide_bounds(binned_counts)

        median = np.median(binned_counts[binned_counts != 0])

        fig,ax = plt.subplots(nrows=2,ncols=1, figsize=(9,6))

        # plot unordered bins
        x1 = np.linspace(0, 1, binned_counts.size)
        ax[0].scatter(x1, binned_counts)
        ax[0].set_yscale("log", base=2)
        ax[0].plot(x1, np.ones(binned_counts.size)*median, c="red", linewidth=3)
        ax[0].plot(x1, np.ones(binned_counts.size)*upper_bound, c="red", linestyle=":", linewidth=3)
        ax[0].plot(x1, np.ones(binned_counts.size)*lower_bound, c="red", linestyle=":", linewidth=3)
        y_lim = ax[0].get_ylim()
        ax[0].set_ylabel("Read Count")
        ax[0].set_xlabel("Unordered Bins")

        # plot model fit
        emp_probs = rbins / rbins.sum()
        m,b = p # slope intercept in log space

        x2 = np.linspace(0, 1, emp_probs.size)
        y = np.array([m*xi + b for xi in x2])
        y = np.power(2, y)
        y /= y.sum()
        ax[1].scatter(x2, emp_probs, c="C1")
        ax[1].plot(x2, y, c="black", linewidth=3)
        ax[1].set_yscale("log", base=2)
        ax[1].set_ylabel("Log2 Probability")
        ax[1].set_xlabel("Reordered Bins")
        ax[1].set_title("\nslope={:.3f} (Reads = {})".format(m, int(estimate.nreads)))

        plt.tight_layout()
        plt.savefig(os.path.join(plot_folder, estimate.sample_id + "-" + estimate.genome_id + ".pdf"))
        plt.close()
        nplots += 1



def estimate_ptrs_coptr_contig(assembly_genome_ids, grouped_coverage_map_folder, min_reads, min_samples, threads, plot_folder=None):
    """Estimate Peak-to-Trough ratios across samples.

    Parameters
    ----------
        assembly_genome_ids : set[str]
            A list of genome ids to estimate PTRs.
        grouped_coverage_map_folder : str
            Path to folder containing coverage maps grouped by genome.
        min_reads : float
            Minumum read count after filtering to estimate PTR
        min_samples : float
            Minimum number of passing samples required to estimate a PTR
        threads : int
            Number of threads for parallel computation
        plot_folder : str
            If specified, plots of fitted model will be saved here

    Returns
    -------
        coptr_contig_estimates : dict[str -> list[CoPTContigEstimate]]
            A dictionary with key reference genome id and value
            a list of CoPTRContigEstimate for that reference genome.
    """
    print_info("CoPTRContig", "checking reference genomes")
    coptr_contig_estimates = {}
    coptr_contig = CoPTRContig(min_reads, min_samples)

    total_reads = 0
    if threads == 1:
        for ref_id in sorted(assembly_genome_ids):
            coverage_maps = []
            with open(os.path.join(grouped_coverage_map_folder, ref_id + ".cm.pkl"), "rb") as f:
                try:
                    while True:
                        coverage_maps.append(pkl.load(f))
                except EOFError:
                    pass

            if plot_folder is not None:
                estimates, parameters, reordered_bins = coptr_contig.estimate_ptrs(coverage_maps, return_bins=True)
                plot_fit(estimates, parameters, coverage_maps, reordered_bins, plot_folder)
            else:
                estimates = coptr_contig.estimate_ptrs(coverage_maps)
            coptr_contig_estimates[ref_id] = estimates

    else:
        # workers for multiprocessing
        pool = mp.Pool(threads)
        coverage_maps = {}

        for i,ref_id in enumerate(sorted(assembly_genome_ids)):
            coverage_maps[ref_id] = []
            with open(os.path.join(grouped_coverage_map_folder, ref_id + ".cm.pkl"), "rb") as f:
                try:
                    while True:
                        coverage_maps[ref_id].append(pkl.load(f))
                except EOFError:
                    pass

            if (i % threads == 0 and i > 0) or i == len(assembly_genome_ids) - 1:
                flat_coverage_maps = [(ref_id, coverage_maps[ref_id], plot_folder) for ref_id in coverage_maps]
                flat_results = pool.map(coptr_contig._parallel_helper, flat_coverage_maps)
                coptr_contig_estimates.update(dict(flat_results))

                coverage_maps = {}

    return coptr_contig_estimates