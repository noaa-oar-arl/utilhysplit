import numpy as np


def example_bootstrap():
    data = np.random.poisson(10, size=100)

    # Number of bootstrap resamples
    B = 10000
    boot_means = np.array([
        np.mean(np.random.choice(data, size=len(data), replace=True))
        for _ in range(B)
    ])

    # 95% confidence interval
    ci_lower, ci_upper = np.percentile(boot_means, [2.5, 97.5])
    print(f"Bootstrap 95% CI for the mean: ({ci_lower:.2f}, {ci_upper:.2f})")

