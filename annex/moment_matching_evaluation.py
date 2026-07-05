import numpy as np
from dataclasses import dataclass
from numpy.typing import NDArray
from scipy.optimize import newton
from scipy.special import i0e, i1e
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import joblib
from joblib import Parallel, delayed
from pqdm.processes import pqdm
import os

rng = np.random.default_rng(seed=1234)


@dataclass
class NormalDist:
    mu: NDArray
    a: NDArray


@dataclass
class GVMDist:
    mu: NDArray
    a: NDArray
    alpha: float
    beta: NDArray
    kappa: float


# First coordinate is assumed to be angular, wrapped to [-pi, pi)
# Returns an array of size (num, dim)
def sample_normal_wrapped(dist: NormalDist, num: int) -> NDArray:
    dim = dist.mu.shape[0]

    z = rng.normal(0.0, 1.0, size=(num, dim))
    x = z @ dist.a.T + dist.mu
    # Wrap
    x[:, 0] = np.mod(x[:, 0] + np.pi, 2.0 * np.pi) - np.pi
    return x


def sample_gvm(dist: GVMDist, num: int) -> NDArray:
    dim = dist.mu.shape[0]

    z_eucl = rng.normal(0.0, 1.0, size=(num, dim))
    x_eucl = z_eucl @ dist.a.T + dist.mu
    thetas = dist.alpha + z_eucl @ dist.beta
    vm = rng.vonmises(thetas, dist.kappa, size=num)
    # Returned samples are in [-theta + pi, theta + pi), i.e. centered around the mean,
    # we displace them so they live in[-pi, pi) centered around 0
    vm = np.mod(vm + np.pi, 2.0 * np.pi) - np.pi

    return np.column_stack((vm, x_eucl))


def moment_match(dist: NormalDist) -> GVMDist:
    # P = P_e
    p = dist.a @ dist.a.T

    # A can be obtained by cholesky afterwards
    a = np.linalg.cholesky(p[1:, 1:])

    gamma = p[1:, 0]
    sigma2 = p[0, 0]

    alpha = dist.mu[0]
    mu = dist.mu[1:]

    # beta = A^(-1) gamma, or in numpy, transposed
    beta = gamma @ np.linalg.inv(a).T

    betabeta = np.dot(beta, beta)

    v = 0.5 * (betabeta - sigma2)
    v = min(v, -1e-15)
    target = np.exp(v)

    if target < 0.999:
        kappa = newton(
            func=lambda k: i1e(k) / i0e(k) - target,
            x0=0.0,
        )
    else:
        # Asymptotic approximation, mostly equal after kappa ~ 500, i.e. target > 0.999
        kappa = 0.5 / (1.0 - target)

    return GVMDist(
        mu=mu,
        a=a,
        alpha=alpha,
        beta=beta,
        kappa=kappa
    )

def train_clasiffier(samples1: NDArray, samples2: NDArray) -> float:
    # Flatten the two arrays, adding the classification
    X = np.vstack((samples1, samples2))
    y = np.concatenate((np.zeros(samples1.shape[0]), np.ones(samples2.shape[0])))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1234)
    rf = RandomForestClassifier(random_state=1234)
    rf.fit(X_train, y_train)

    return rf.score(X_test, y_test)

def controlled_cholesky(dim, a00, ai0_norm, max_other_variances=3.0):
    a = np.zeros((dim, dim))
    a[0, 0] = a00

    couple_vector = np.random.randn(dim - 1)
    couple_vector /= np.linalg.norm(couple_vector)
    couple_vector *= ai0_norm

    a[1:, 0] = couple_vector

    random_tril = np.tril(np.random.randn(dim - 1, dim - 1))
    np.fill_diagonal(random_tril, np.abs(np.diagonal(random_tril)))
    row_norms = np.linalg.norm(random_tril, axis=1, keepdims=True)
    normalized_tril = random_tril / row_norms

    target_vars = np.random.uniform(1.0, max_other_variances, size=dim - 1)
    target_vars = np.maximum(target_vars, couple_vector**2 + 0.1)
    remaining_std = np.sqrt(target_vars - couple_vector**2)
    a[1:, 1:] = normalized_tril * remaining_std[:, None]

    return a

def evaluate_grid(a00, ai0, n_samples, n_average, dim_euclid):
    out = {}

    out["a00"] = a00
    out["ai0"] = ai0
    other_variances = 10.0
    out["other_variances"] = other_variances

    tries = []

    for ntry in range(0, n_average):
        out_try = {}
        
        a = controlled_cholesky(dim_euclid + 1, a00, ai0, other_variances)
        out_try["a"] = a

        mu = np.random.normal(0.0, 10.0, dim_euclid + 1)
        out_try["mu"] = mu

        norm = NormalDist(mu, a)
        gvm = moment_match(norm)

        out_try["gvm_a"] = gvm.a
        out_try["gvm_alpha"] = gvm.alpha
        out_try["gvm_beta"] = gvm.beta
        out_try["gvm_kappa"] = gvm.kappa
        out_try["gvm_mu"] = gvm.mu

        norm_samples = sample_normal_wrapped(norm, n_samples)
        out_try["norm_samples"] = norm_samples

        gvm_samples = sample_gvm(gvm, n_samples)
        out_try["gvm_samples"] = gvm_samples

        out_try["score"] = train_clasiffier(norm_samples, gvm_samples)

        tries.append(out_try)
    
    out["tries"] = tries
    return out


if __name__ == "__main__":
    # We sweep A[0,0] (sigma) and norm(A[i,0]) (coupling between linear and angular dimension) in a grid
    sweep_min = (1e-6, 1e-6)
    sweep_max = (6.0, 6.0)

    # Number of samples per training
    n_samples = 5000
    # Number of euclidean dimensions on top of the angular
    dim_euclid = 3

    # Number of points in the grid in each direction
    n_points00 = 10
    n_pointsi0 = 10

    # Number of averages per grid square
    n_average = 5

    a00vals = np.linspace(sweep_min[0], sweep_max[0], n_points00)
    ai0vals = np.linspace(sweep_min[1], sweep_max[1], n_pointsi0)

    tasks = [
        (a00, a10, n_samples, n_average, dim_euclid)
        for a10 in ai0vals
        for a00 in a00vals
    ]

    results = pqdm(
        tasks, 
        evaluate_grid,  # ty:ignore[invalid-argument-type]
        n_jobs=os.cpu_count() or 1, 
        argument_type='args',
        desc="Sweeping Parameter Grid"
    )

    serialize = {
        "a00vals": a00vals,
        "a10vals": ai0vals,
        "results": results,
    }

    joblib.dump(serialize, "results.pkl")



