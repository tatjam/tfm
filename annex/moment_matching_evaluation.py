import numpy as np
from dataclasses import dataclass
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.special import i0e, i1e

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

    # Solve the transcendental equation for kappa, we could feed the derivative
    # but the secant method just works
    target = np.exp(0.5 * (betabeta - sigma2))
    kappa = newton(
        func=lambda k: (i1e(k) / i0e(k)) - target,
        x0=1.0,
    )

    return GVMDist(
        mu=mu,
        a=a,
        alpha=alpha,
        beta=beta,
        kappa=kappa
    )


dnormal = NormalDist(
    mu=np.array([1.0, 3.0]),
    a=np.array([[0.1, 0.0], [0.0, 1.0]]),
)

dgvm = moment_match(dnormal)

samples = sample_normal_wrapped(dnormal, 100000)
samples2 = sample_gvm(dgvm, 100000)

plt.axis('equal')
plt.scatter(samples[:, 0], samples[:, 1], alpha=0.2, s=4, color='purple', label='Samples')
plt.scatter(samples2[:, 0], samples2[:, 1], alpha=0.2, s=4, color='blue', label='Samples')
plt.show()
