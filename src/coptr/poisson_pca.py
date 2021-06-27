"""
poisson_pca.py
======================
Dimensionality reduction with Poisson loss function.
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
import scipy.optimize


logger = logging.getLogger(__name__)


class PoissonPCA:
    """Principal Component Analysis with Poisson observations."""

    def pca(self, X, k, tol=1e-3):
        """Run the PCA. For a data matrix X (sample by observations), finds two matrices W and V
        such that:

            X_{ij} ~ Poisson(exp( [WV]_{ij} ))

        Note that missing entries are possible. Missing entries should be denoted by np.nan.

        Parameters
        ----------
            X : np.array
                A sample by observation data matrix with observed counts
            k : int
                Number of latent dimensions. Determines the rank of W and V.
        """
        assert k <= X.shape[1], "k must be smaller than the number of columns of X"
        assert k > 0 and k == int(k), "k must be a positive integer"
        assert np.all(
            X[np.isfinite(X)] >= 0
        ), "X must have nonnegative entries or be missing"
        X = X.astype(float)

        # use PCA to initialize matrices W,V
        # initializing missing values to mean
        Xt = np.copy(X)
        Xt[Xt == 0] = 0.1
        for i, row in enumerate(Xt):
            if np.isnan(row).sum() > 0:
                Xt[i, np.isnan(row)] = np.mean(row[np.isfinite(row)])
        u, s, vh = np.linalg.svd(np.log2(Xt))
        u = u[:, :k]
        s = s[:k]
        vh = vh[:k, :]
        W = u * s
        V = vh

        it = 0
        prv_log_lk = -np.inf
        log_lk = self.log_likelihood(X, W, V)
        while np.abs(prv_log_lk - log_lk) > tol:

            logger.debug("Iteration: %d; Log likelihood: %.3G", it, log_lk)

            prv_log_lk = log_lk
            W = self.update_W(X, W, V)
            V = self.update_V(X, W, V)
            log_lk = self.log_likelihood(X, W, V)

            it += 1

            assert (
                prv_log_lk <= log_lk or np.abs(prv_log_lk - log_lk) < 1e-6
            ), "log likelihood must be increasing ({} > {})".format(prv_log_lk, log_lk)
        return W, V

    def update_W(self, X, W, V):
        """Optimize the log-likelihood with respect to W.

        Parameters
        ----------
            X : np.array
                The data matrix
            W : np.array
                The current estimate of W
            V : np.array
                The current estimate of V
        """
        loss_fn = lambda Wx: -self.log_likelihood(
            X, Wx.reshape(X.shape[0], V.shape[0]), V
        )
        grad_fn = lambda Wx: -self.grad_W(
            X, Wx.reshape(X.shape[0], V.shape[0]), V
        ).flatten()
        # these are handled by the minimizer
        np.seterr(over="ignore")
        res = scipy.optimize.minimize(
            loss_fn, np.copy(W).flatten(), jac=grad_fn, method="CG"
        )
        np.seterr(over="warn")
        return res.x.reshape(X.shape[0], V.shape[0])

    def update_V(self, X, W, V):
        """Optimize the log-likelihood with respect to V.

        Parameters
        ----------
            X : np.array
                The data matrix
            W : np.array
                The current estimate of W
            V : np.array
                The current estimate of V
        """
        return self.update_W(X.T, V.T, W.T).T

    def log_likelihood(self, X, W, V):
        """Compute the log likelihood (minus a constant).

        Parameters
        ----------
            X : np.array
                An n x m data matrix
            W : np.array
                An n x k matrix
            V : np.array
                A k x m matrix
        """
        param = np.exp(W.dot(V))
        return (X * W.dot(V) - param)[np.isfinite(X)].sum()

    def grad_W(self, X, W, V):
        """Compute the gradient of the log likelihood wrt W.

        Parameters
        ----------
            X : np.array
                An n x m data matrix
            W : np.array
                An n x k matrix
            V : np.array
                A k x m matrix
        """
        Xtmp = np.copy(X)
        # zero out rows so they do not contribute to inner product
        Xtmp[np.isnan(X)] = 0
        grad = Xtmp.dot(V.T) - np.exp(W.dot(V)).dot(V.T)
        return grad


if __name__ == "__main__":
    np.random.seed(101441)

    ppca = PoissonPCA()

    # generate 100 random test sets
    for i in range(100):
        print("running test set", i + 1, "of", 100)
        n = np.random.randint(10, 100)
        m = np.random.randint(5, n)
        k = np.random.randint(1, np.min((m, 10)))
        print("\tn {} m {} k {}".format(n, m, k))

        # generate test data
        W = np.random.normal(loc=0, scale=1, size=(n, k))
        V = np.random.normal(loc=0, scale=1, size=(k, m))

        X = np.random.poisson(np.exp(W.dot(V))).astype(float)

        coords = [(i, j) for i in range(n) for j in range(m)]
        nmissing = int(0.25 * X.size)
        for l in range(nmissing):
            coord = coords[np.random.randint(len(coords))]
            X[coord] = np.nan

        X_finite = X[np.isfinite(X)]
        print("\tmin {}".format(X_finite.min()))
        print("\tmax {}".format(X_finite.max()))
        print("\tmean {}".format(X_finite.mean()))
        print("\tstd {}".format(X_finite.std()))

        West, Vest = ppca.pca(X, k)

    print("passed!")
