"""
coptr_contig.py
======================
Estimate peak-to-trough ratios from assemblies.
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
import multiprocessing as mp
import os.path
import pickle as pkl
import struct
import sys

import numpy as np
import scipy.special
import scipy.stats

from .poisson_pca import PoissonPCA


logger = logging.getLogger(__name__)


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
        return "CoPTRContigEstimate(bam_file={}, genome_id={}, sample_id={}, estimate={:.3f}, nreads={}, cov_frac={})".format(
            self.bam_file,
            self.genome_id,
            self.sample_id,
            self.estimate,
            self.nreads,
            self.cov_frac,
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

    def compute_genomewide_bounds(self, binned_counts, crit_region=0.05):
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
        m = norm.ppf(1 - crit_region)

        lower_bound = np.power(2, median - m * std)
        upper_bound = np.power(2, median + m * std)

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
        contig_ids = coverage_maps[0].contig_ids
        ref_genome = coverage_maps[0].genome_id

        for cm in coverage_maps:

            binned_reads = []
            for contig_id in sorted(contig_ids):

                if contig_id not in cm.contig_ids:
                    logger.critical(
                        "Missing contig %s from %s in %s",
                        contig_id,
                        ref_genome,
                        cm.sample_id,
                    )
                    sys.exit(1)

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
        x[:, 0] = rng

        sol, residuals, rank, s = np.linalg.lstsq(x, y, rcond=None)
        m, b = sol

        pred = m * rng + b
        residuals = y - pred

        return m, b, residuals

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
        # change to base e
        loge_ptr = log2_ptr / np.log2(np.e)
        log_probs = np.array([loge_ptr * xi for xi in x])
        log_norm = scipy.special.logsumexp(log_probs)
        log_probs -= log_norm

        # x = np.linspace(0, 1, binned_counts.size)
        # # compute normalized probabilities
        # log_probs = np.array([log2_ptr*xi for xi in x])
        # log_norm = np.log2(np.power(2, log_probs).sum())
        # log_probs -= log_norm

        log_lk = (binned_counts * log_probs).sum()

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

        log2_ptr, fun, it, funcalls = scipy.optimize.brent(f, full_output=True)
        log_lk = -fun

        x = np.linspace(0, 1, binned_counts.size)
        # compute normalized probabilities
        loge_ptr = log2_ptr / np.log2(np.e)
        log_probs = np.array([loge_ptr * xi for xi in x])
        log_norm = scipy.special.logsumexp(log_probs)
        log_probs -= log_norm
        probs = np.exp(log_probs)
        intercept = probs[0]

        # the orientation of the bins is arbitrary,
        # so let's orient with a positive slope
        flip = False
        if log2_ptr < 0:
            intercept = log2_ptr + intercept
            log2_ptr = np.abs(log2_ptr)
            flip = True

        return log2_ptr, intercept, log_lk, flip

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
                estimate = CoPTRContigEstimate(
                    cm.bam_file,
                    cm.genome_id,
                    cm.sample_id,
                    np.nan,
                    cm.reads,
                    cm.frac_nonzero,
                )
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
        for i, col in enumerate(A.T):

            lower_bound, upper_bound = self.compute_genomewide_bounds(
                col, crit_region=0.05
            )
            col[np.logical_or(col > upper_bound, col < lower_bound)] = np.nan

            frac_filtered = np.isnan(col).sum() / col.size
            reads = col[np.isfinite(col)].sum()
            if reads > self.min_reads and frac_filtered < 0.25:
                A_filtered.append(col)
                bam_files.append(passing_coverage_maps[i].bam_file)
                genome_ids.append(passing_coverage_maps[i].genome_id)
                sample_ids.append(passing_coverage_maps[i].sample_id)
                tmp.append(passing_coverage_maps[i])
            else:
                bam_file = passing_coverage_maps[i].bam_file
                genome_id = passing_coverage_maps[i].genome_id
                sample_id = passing_coverage_maps[i].sample_id
                frac_nonzero = passing_coverage_maps[i].frac_nonzero
                m = np.nan
                estimate = CoPTRContigEstimate(
                    bam_file, genome_id, sample_id, np.nan, reads, frac_nonzero
                )
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])

        if len(A_filtered) < self.min_samples:
            for cm in tmp:
                bam_file = cm.bam_file
                genome_id = cm.genome_id
                estimate = CoPTRContigEstimate(
                    cm.bam_file,
                    cm.genome_id,
                    cm.sample_id,
                    np.nan,
                    cm.reads,
                    cm.frac_nonzero,
                )
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])
            if return_bins:
                return estimates, parameters, binned_counts
            else:
                return estimates
        else:
            tmp = None

        logger.info("Running %s.", genome_id)

        A = np.array(A_filtered).T

        # filter out rows where too many bins are missing
        keep_rows = np.sum(np.isnan(A), axis=1) < 0.025 * A.shape[1]
        A = A[keep_rows, :]

        ppca_failed = False
        try:
            if keep_rows.sum() >= 0.5 * keep_rows.size:
                poisson_pca = PoissonPCA()
                W, V = poisson_pca.pca(A, k=1)
        except:
            ppca_failed = True
            logger.warn("PoissonPCA failed for assembly " + genome_ids[0])

        if keep_rows.sum() < 0.5 * keep_rows.size or ppca_failed:
            for j, col in enumerate(A.T):
                reads = col[np.isfinite(col)].sum()
                frac_nonzero = (col[np.isfinite(col)] > 0).sum() / col[
                    np.isfinite(col)
                ].size
                estimate = CoPTRContigEstimate(
                    bam_files[j],
                    genome_ids[j],
                    sample_ids[j],
                    np.nan,
                    reads,
                    frac_nonzero,
                )
                estimates.append(estimate)
                parameters.append((np.nan, np.nan))
                binned_counts.append([np.nan])

            if return_bins:
                return estimates, parameters, binned_counts
            else:
                return estimates

        sorted_idx = np.argsort(W.flatten())
        A = A[sorted_idx, :]

        for j, col in enumerate(A.T):
            col = col[np.isfinite(col)]
            reads = col.sum()
            lb = int(0.05 * col.size)
            ub = int(0.95 * col.size)
            col = col[lb:ub]

            m, b, log_lk, flip = self.estimate_ptr_maximum_likelihood(col)
            # m,b,res = self.estimate_slope_intercept(col)
            # flip = m < 0

            # the sign of the slope depends on the orientation,
            # so let's look at positive slopes only
            m = np.abs(m)
            frac_nonzero = (col > 0).sum() / col.size
            estimate = CoPTRContigEstimate(
                bam_files[j], genome_ids[j], sample_ids[j], m, reads, frac_nonzero
            )
            estimates.append(estimate)

            parameters.append((m, b))

            bc = col
            if flip:
                bc = np.flip(bc)
            binned_counts.append(bc)

        logger.info("Finished %s.", genome_id)

        if return_bins:
            return estimates, parameters, binned_counts
        else:
            return estimates

    def _parallel_helper(self, x):
        ref_genome = x[0]
        coverage_maps = x[1]
        plot_folder = x[2]

        if plot_folder is not None:
            estimates, parameters, reordered_bins = self.estimate_ptrs(
                coverage_maps, return_bins=True
            )
            plot_fit(estimates, parameters, coverage_maps, reordered_bins, plot_folder)
        else:
            estimates = self.estimate_ptrs(coverage_maps)

        return (ref_genome, estimates)


def plot_fit(estimates, parameters, coverage_maps, reordered_bins, plot_folder):
    import matplotlib.pyplot as plt

    coptr_contig = CoPTRContig(min_reads=5000, min_samples=5)

    # make sure esimates and coverage maps are in the same order
    order = dict([(est.sample_id, i) for i, est in enumerate(estimates)])
    tmp = [[] for cm in coverage_maps]
    for cm in coverage_maps:
        idx = order[cm.sample_id]
        tmp[idx] = cm
    coverage_maps = tmp

    nplots = 0
    for estimate, p, cm, rbins in zip(
        estimates, parameters, coverage_maps, reordered_bins
    ):
        if np.isnan(estimate.estimate):
            continue

        contig_ids = cm.contig_ids
        lengths = []
        binned_counts = []
        for cid in contig_ids:
            length = cm.get_length(cid)
            if length < 11000:
                continue

            binned_counts += cm.bin_reads(cid).tolist()

        binned_counts = np.array(binned_counts)
        lower_bound, upper_bound = coptr_contig.compute_genomewide_bounds(binned_counts)

        median = np.median(binned_counts[binned_counts != 0])

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))

        # plot unordered bins
        x1 = np.linspace(0, 1, binned_counts.size)
        ax[0].scatter(x1, binned_counts)
        ax[0].set_yscale("log", base=2)
        ax[0].plot(x1, np.ones(binned_counts.size) * median, c="red", linewidth=3)
        ax[0].plot(
            x1,
            np.ones(binned_counts.size) * upper_bound,
            c="red",
            linestyle=":",
            linewidth=3,
        )
        ax[0].plot(
            x1,
            np.ones(binned_counts.size) * lower_bound,
            c="red",
            linestyle=":",
            linewidth=3,
        )
        y_lim = ax[0].get_ylim()
        ax[0].set_ylabel("Read Count")
        ax[0].set_xlabel("Unordered Bins")

        # plot model fit
        emp_probs = rbins / rbins.sum()
        m, b = p  # slope intercept in log space

        x2 = np.linspace(0, 1, emp_probs.size)
        y = np.array([m * xi + b for xi in x2])
        y = np.power(2, y)
        y /= y.sum()
        ax[1].scatter(x2, emp_probs, c="C1")
        ax[1].plot(x2, y, c="black", linewidth=3)
        ax[1].set_yscale("log", base=2)
        ax[1].set_ylabel("Density")
        ax[1].set_xlabel("Reordered Bins")
        ax[1].set_title(
            "\nlog2(PTR)={:.3f} (Reads = {})".format(m, int(estimate.nreads))
        )

        plt.tight_layout()
        plt.savefig(
            os.path.join(
                plot_folder, estimate.sample_id + "-" + estimate.genome_id + ".pdf"
            )
        )
        plt.close()
        nplots += 1


def load_coverage_maps(file_handle, ref_id):
    reached_end = False
    coverage_maps = []
    try:
        while True:
            cm = pkl.load(file_handle)
            if cm.passed_qc():
                coverage_maps.append(cm)

    except EOFError:
        pass
    return coverage_maps


def estimate_ptrs_coptr_contig(
    assembly_genome_ids,
    grouped_coverage_map_folder,
    min_reads,
    min_samples,
    plot_folder=None,
):
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
    logger.info("Checking reference genomes.")
    coptr_contig_estimates = {}
    coptr_contig = CoPTRContig(min_reads, min_samples)

    for ref_id in sorted(assembly_genome_ids):
        with open(
            os.path.join(grouped_coverage_map_folder, ref_id + ".cm.pkl"), "rb"
        ) as f:
            coverage_maps = load_coverage_maps(f, ref_id)

        if plot_folder is not None:
            estimates, parameters, reordered_bins = coptr_contig.estimate_ptrs(
                coverage_maps, return_bins=True
            )
            plot_fit(
                estimates, parameters, coverage_maps, reordered_bins, plot_folder
            )
        else:
            estimates = coptr_contig.estimate_ptrs(coverage_maps)

        coptr_contig_estimates[ref_id] = estimates

    return coptr_contig_estimates
