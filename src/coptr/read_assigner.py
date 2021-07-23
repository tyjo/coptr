"""
read_assigner.py
======================
Assign multi-mapped reads to a single genome.
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

import numpy as np
import scipy.misc
import scipy.special
import scipy.stats
from scipy.sparse import csr_matrix
from scipy.special import digamma


logger = logging.getLogger(__name__)


class ReadAssigner(object):
    """Assign multi-mapped reads using a mixture model.

    Parameters
    ----------
        X : np.array
            A read by genome data matrix (0-1). X[i,j] = 1
            if read i maps to genome j
    """

    def __init__(self, X, prior_counts):
        self.X = csr_matrix(X)
        self.nreads = X.shape[0]
        self.ngenomes = X.shape[1]
        assert np.all(
            self.X.sum(axis=1) > 0
        ), "each read must have at least one valid mapping"

        # set alpha based on unambiguously assigned reads
        self.alpha = np.ones(self.ngenomes) + prior_counts

        # variational parameters
        self.phi = self.X.multiply(1 / self.X.sum(axis=1))
        self.eta = np.copy(self.alpha)

    def compute_elbo(self):
        """Compute the variational objective function."""
        elbo = 0

        # E_q[log p(pi)]
        elbo += ((self.alpha - 1) * (digamma(self.eta) - digamma(self.eta.sum()))).sum()

        # E_q[log p(X | Z)]
        elbo += (self.phi.multiply(self.X)).sum()

        # E_q[log p(Z)]
        elbo += (self.phi.dot(digamma(self.eta) - digamma(self.eta.sum()))).sum()

        # entropy of dirichlet
        elbo += scipy.stats.dirichlet(self.eta).entropy()

        # entropy of categorical
        tmp = self.phi.data * np.log(self.phi.data)
        elbo += -tmp.sum()

        return elbo

    def update_phi(self):
        """Optimize phi with fixed eta."""
        log_phi = self.X.multiply(digamma(self.eta))
        phi = csr_matrix(
            (np.exp(log_phi.data), (log_phi.row, log_phi.col)), shape=log_phi.shape
        )
        self.phi = phi.multiply(1 / phi.sum(axis=1))

    def update_eta(self):
        """Optimize eta with fixed phi."""
        self.eta = np.array(self.alpha + self.phi.sum(axis=0))[0]

    def run_vi(self, tol=1e-3):
        """Optimize variational parameters with respect to elbo."""
        prv_elbo = -np.inf
        elbo = self.compute_elbo()

        it = 0
        while np.abs(prv_elbo - elbo) > tol:
            logger.debug("Iteration: %d; elbo: %.3G", it, elbo)
            self.update_eta()
            self.update_phi()
            prv_elbo = elbo
            elbo = self.compute_elbo()
            it += 1

            assert (
                elbo >= prv_elbo - 1e-3
            ), "elbo must be strictly increasing ({} > {})".format(elbo, prv_elbo)

        logger.debug("Iteration: %d; elbo: %.3G", it, elbo)

    def assign_reads(self):
        """Compute read assignments.

        Return
        ------
            assignments : np.array
                A 1D array. Each entry gives a genome id that is
                the genome assignment for each read.
        """
        self.run_vi()
        self.phi = csr_matrix(self.phi)
        assignments = []
        for i in range(self.nreads):
            assignments.append(np.argmax(self.phi[i].toarray()))
        assignments = np.array(assignments)
        return assignments


if __name__ == "__main__":
    np.random.seed(103254)

    # test a large number of cases
    for i in range(100):
        print("running test set", i + 1)
        nreads = np.random.randint(10000, 100000)
        ngenomes = np.random.randint(10, 50)
        Z = np.zeros((nreads, ngenomes))

        print("\treads", nreads)
        print("\tgenomes", ngenomes)

        for n in range(nreads):
            true_assignment = np.random.randint(5)
            Z[n, true_assignment] = 1

        X = np.copy(Z)
        multimapped_reads = int(0.1 * nreads)
        for n in range(multimapped_reads):
            true_assignment = np.argwhere(Z[n] == 1).flatten()[0]
            additional_mappings = np.max((1, np.random.poisson(1)))
            for i in range(additional_mappings):
                mapping_id = np.random.randint(ngenomes)
                X[n, mapping_id] = 1
        X = X[:multimapped_reads, :]

        prior = Z[
            multimapped_reads:,
        ].sum(axis=0)
        read_assigner = ReadAssigner(X[:multimapped_reads], prior)
        assignments = read_assigner.assign_reads()
        correct_assignments = 0
        for i, assignment in enumerate(assignments):
            if Z[i, assignment] == 1:
                correct_assignments += 1
        print(
            "\t{} of reads correctly mapped".format(
                correct_assignments / len(assignments)
            )
        )
    print("done!")
