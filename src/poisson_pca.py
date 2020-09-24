import numpy as np
import scipy.optimize
import sys

class PoissonPCA:
    """Principal Component Analysis with Poisson observations.
    """


    def pca(self, X, k, tol=1E-3, verbose=False):
        """Run the PCA. For a data matrix X (sample by observations), finds two matrices W and V
        such that:

            X_{ij} ~ Poisson(exp( [WV]_{ij} ))

        Parameters
        ----------
            X : np.array
                A sample by observation data matrix with observed counts
            k : int
                Number of latent dimensions. Determines the rank of W and V.
        """
        assert k <= X.shape[1], "k must be smaller than the number of columns of X"
        assert k > 0 and k == int(k), "k must be a positive integer"
        assert np.all(X >= 0), "X must have nonnegative entries"
        X = X.astype(float)

        # use PCA to initialize matrices W,V
        Xt = np.copy(X)
        Xt[Xt == 0] = 0.1
        u,s,vh = np.linalg.svd(np.log2(Xt))
        u = u[:,:k]
        s = s[:k]
        vh = vh[:k,:]
        W = u*s
        V = vh

        it = 0
        prv_log_lk = -np.inf
        log_lk = self.log_likelihood(X, W, V)
        while np.abs(prv_log_lk - log_lk) > tol:

            if verbose:
                print("iteration:", it, "log likelihood:", log_lk, file=sys.stderr)

            prv_log_lk = log_lk
            W = self.update_W(X, W, V)
            V = self.update_V(X, W, V)
            log_lk = self.log_likelihood(X, W, V)

            it += 1

            assert prv_log_lk <= log_lk or np.abs(prv_log_lk - log_lk) < 1E-6, "log likelihood must be increasing ({} > {})".format(prv_log_lk, log_lk)
        return W,V


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
        loss_fn = lambda Wx : -self.log_likelihood(X, Wx.reshape(X.shape[0], V.shape[0]), V)
        grad_fn = lambda Wx : -self.grad_W(X, Wx.reshape(X.shape[0], V.shape[0]), V).flatten()
        # these are handled by the minimizer
        np.seterr(overflow="ignore")
        res = scipy.optimize.minimize(loss_fn, np.copy(W).flatten(), jac=grad_fn, method="CG")
        np.seterr(overflow="warn")
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


    def update_WV(self, X, W, V):
        def loss_fn(WVT_flat):
            W_est = WVT_flat[:W.size].reshape(W.shape)
            V_est = WVT_flat[W.size:].reshape(V.T.shape).T
            return -self.log_likelihood(X, W_est, V_est)

        def grad_fn(WVT_flat):
            W_est = WVT_flat[:W.size].reshape(W.shape)
            V_est = WVT_flat[W.size:].reshape(V.T.shape).T
            grad_W = self.grad_W(X, W_est, V_est)
            grad_V = self.grad_V(X, W_est, V_est)
            return np.vstack((grad_W, grad_V.T)).flatten()

        WV = np.vstack((W, V.T))
        res = scipy.optimize.minimize(loss_fn, WV.flatten(), jac=grad_fn, method="CG")
        #print(res.fun)
        #quit()
        WV = res.x.reshape(WV.shape)
        W = WV[:W.shape[0]]
        V = WV[W.shape[0]:].T
        return W,V



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
        return (X*W.dot(V) - param).sum()


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
        return X.dot(V.T) - np.exp(W.dot(V)).dot(V.T)
        


    def grad_V(self, X, W, V):
        """Compute the gradient of the log likelihood wrt V.

        Parameters
        ----------
            X : np.array
                An n x m data matrix
            W : np.array
                An n x k matrix
            V : np.array
                A k x m matrix
        """
        return self.grad_W(X.T, V.T, W.T).T


if __name__ == "__main__":
    np.random.seed(101441)

    ppca = PoissonPCA()

    # generate 100 random test sets
    for i in range(100):
        print("running test set", i+1, "of", 100)
        n = np.random.randint(10, 100)
        m = np.random.randint(5, n)
        k = np.random.randint(1, np.min((m, 10)))
        print("\tn {} m {} k {}".format(n,m,k))

        # generate test data
        W = np.random.normal(loc=0, scale=1, size=(n, k))
        V = np.random.normal(loc=0, scale=1, size=(k, m))

        X = np.random.poisson(np.exp(W.dot(V)))
        print("\tmin {}".format(X.min()))
        print("\tmax {}".format(X.max()))
        print("\tmean {}".format(X.mean()))
        print("\tstd {}".format(X.std()))

        West, Vest = ppca.factor(X, k)

    print("passed!")
