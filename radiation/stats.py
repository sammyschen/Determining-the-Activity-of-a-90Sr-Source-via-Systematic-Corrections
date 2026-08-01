"""Small weighted-statistics helpers shared across the analysis scripts."""
import numpy as np


def weighted_mean(values, sigmas):
    """Weighted mean using 1/sigma^2 weights. Returns (mean, sigma_of_mean)."""
    values = np.asarray(values, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    w = 1.0 / sigmas**2
    mean = np.sum(w * values) / np.sum(w)
    sigma = np.sqrt(1.0 / np.sum(w))
    return mean, sigma


def weighted_linear_fit(x, y, sigma_y):
    """Weighted least squares fit y = m*x + c.

    Returns (m, sigma_m, c, sigma_c).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = 1.0 / np.asarray(sigma_y, dtype=float) ** 2
    W = np.diag(w)
    X = np.vstack([x, np.ones_like(x)]).T
    XTWX = X.T @ W @ X
    XTWy = X.T @ W @ y
    m, c = np.linalg.solve(XTWX, XTWy)
    cov = np.linalg.inv(XTWX)
    return m, np.sqrt(cov[0, 0]), c, np.sqrt(cov[1, 1])
