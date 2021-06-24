import pickle as pkl
import unittest

import numpy as np
import scipy.stats

from src.coptr.bam_processor import CoverageMapContig
from src.coptr.coptr_contig import CoPTRContig
from src.coptr.coptr_ref import CoPTRRef


class TestCoPTR(unittest.TestCase):
    def test_coptr(self):
        print("running test simulations for CoPTRRef and CoPTRContig...")
        np.random.seed(22721)

        simulator = Simulator()
        coptr_ref = CoPTRRef(min_reads=5000, min_cov=0.75)
        coptr_contig = CoPTRContig(min_reads=5000, min_samples=5)
        coptr_ref_results = []
        coptr_contig_results = []

        # simulate 100 samples
        sim_log2_ptrs = []
        coverage_map_contig = []
        ori_pos = np.random.random()
        for i in range(20):
            print("\tsimulation", i + 1, "of", 20)
            read_positions, ptr, ori_pos, ter_pos = simulator.simulate_reads(
                20000, ori_pos=ori_pos
            )
            read_positions = read_positions.astype(int)
            sim_log2_ptrs.append(np.log2(ptr))

            est_log2_ptr, ori, ter, best_f, qc_result = coptr_ref.estimate_ptr(
                read_positions, simulator.get_genome_length()
            )
            coptr_ref_results.append(est_log2_ptr)

            coverage_map_contig.append(
                CoverageMapContig(
                    "{}.bam".format(i),
                    "E-coli-simulation",
                    {"e-coli": read_positions},
                    {"e-coli": simulator.get_genome_length()},
                )
            )

            if len(coverage_map_contig) % 10 == 0:
                coptr_contig_estimates = sorted(
                    coptr_contig.estimate_ptrs(coverage_map_contig),
                    key=lambda e: e.bam_file,
                )
                coptr_contig_results += [e.estimate for e in coptr_contig_estimates]
                coverage_map_contig = []
                ori_pos = np.random.random()

        corr_ref = scipy.stats.pearsonr(coptr_ref_results, sim_log2_ptrs)[0]
        corr_contig = scipy.stats.pearsonr(coptr_contig_results, sim_log2_ptrs)[0]

        self.assertTrue(corr_ref > 0.9)
        self.assertTrue(corr_contig > 0.9)


class Simulator:
    """Simulate sequencing reads using the density of reads along an E. coli genome."""

    def __init__(self):
        f = open("test/e-coli-NZ_CP011495.1.pkl", "rb")
        read_density_map, ptr = pkl.load(f)
        f.close()

        self.bin_edges = read_density_map[:, 0]
        self.probs = read_density_map[:, 1]  # left bin edges
        self.min_log2_ptr = 0.05
        self.max_log2_ptr = 1.25

    def get_genome_length(self):
        """Return the number of bases in the genome."""
        return 100 * self.probs.size

    def simulate_reads(self, nreads, ptr=None, ori_pos=None):
        """Simulate nreads along 100bp bins with specified ptr and
        origin position. If ptr or ori_pos are None, they are chosen
        randomly.

        Parameters
        ----------
            nreads : int
                Number of reads to simulate
            ptr : float
                PTR to simulate. If None, a PTR is chosen at random.
            ori_pos : float
                In the interval [0, 1]. Position of the replication origin.
                If None, a position is chosen at random.

        Returns
        -------
            read_positions : np.array
                Coordinates of simulated reads along the genome
            ptr : float
                The simulated PTR
            ori_pos : float
                The simulated replication origin as a fraction along the genome
            ter_pos : float
                The simulated replication terminus position as a fraction along the genome
        """
        if ptr is None:
            ptr = (
                np.random.random() * (self.max_log2_ptr - self.min_log2_ptr)
                + self.min_log2_ptr
            )
            ptr = 2 ** ptr
        else:
            assert np.log2(ptr) > self.min_log2_ptr, "ptr < 2^0.5 may be unreliable"

        if ori_pos is None:
            ori_pos = np.random.random()
        else:
            assert ori_pos >= 0 and ori_pos <= 1, "oriC coordinates must be in [0, 1]"

        ter_pos = (ori_pos + 0.5) % 1
        adj_probs = self.adjust_read_probs(ptr, ori_pos, ter_pos)

        # 100bp bins
        binned_counts = np.random.multinomial(nreads, adj_probs)

        # convert counts to coordinates
        read_positions = []
        for i, c in enumerate(binned_counts):
            size = binned_counts.size
            m = 0.5 * (2 * i + 1) / size
            read_positions += [m for i in range(c)]
        read_positions = self.get_genome_length() * np.array(read_positions)

        return read_positions, ptr, ori_pos, ter_pos  # 100bp bins

    def adjust_read_probs(self, ptr, ori_pos, ter_pos):
        """Scale bin probabilities based on the PTR.

        Parameters
        ----------
            ptr : float
                The PTR to simulate
            ori_pos : float
                In [0,1]. The position of the replication origin
            ter_pos : float
                In [0,1]. The position of the replication terminus
        """
        alpha = np.log2(ptr) / (ori_pos - ter_pos)

        adj_probs = np.zeros(self.probs.size)
        for i, p in enumerate(self.probs):

            # we're approximating an integral over
            # (bin_edges[i], bin_edges[i+1]), so let's
            # take the midpoint of each bin
            if i == self.bin_edges.size - 1:
                m = 0.5 * (self.bin_edges[i] + 1)
                length = 1 - self.bin_edges[i]
            else:
                m = 0.5 * (self.bin_edges[i] + self.bin_edges[i + 1])
                length = self.bin_edges[i + 1] - self.bin_edges[i]

            x1 = np.min([ori_pos, ter_pos])
            x2 = np.max([ori_pos, ter_pos])

            # since we normalize bins at the end,
            # we can compute c1 and c2 + a constant
            # in this case we can use c - log p(x_t)
            if ori_pos < ter_pos:
                c1 = np.log2(ptr)
                c2 = 0
            else:
                c1 = 0
                c2 = np.log2(ptr)

            if m <= x1:
                adj_probs[i] = -alpha * (m - x1) + c1
            elif x1 < m and m < x2:
                adj_probs[i] = alpha * (m - x1) + c1
            else:
                adj_probs[i] = -alpha * (m - x2) + c2
            adj_probs[i] += np.log2(length)

        # normalize
        adj_probs -= np.log2(np.power(2, adj_probs).sum())

        # to log base e
        adj_probs = adj_probs / np.log2(np.exp(1))

        # reweight and normalize
        new_probs = np.zeros(self.probs.size)
        new_probs[self.probs != 0] = (
            np.log(self.probs[self.probs != 0]) + adj_probs[self.probs != 0]
        )
        new_probs[self.probs != 0] -= scipy.special.logsumexp(
            new_probs[self.probs != 0]
        )
        new_probs[self.probs != 0] = np.exp(new_probs[self.probs != 0])

        return new_probs
