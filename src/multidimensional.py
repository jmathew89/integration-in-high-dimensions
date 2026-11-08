import numpy as np
import time
import matplotlib.pyplot as plt

def f(x):
    # x has shape (number of points, dimension)
    return np.sum(x**2, axis=1)

def midpoint_rule(d, m):
    # m midpoint samples along each dimension
    coords = (np.arange(m) + 0.5) / m
    grids = np.meshgrid(*([coords] * d), indexing="ij")
    points = np.stack([g.ravel() for g in grids], axis=1)
    return np.mean(f(points))

def monte_carlo(d, n, rng):
    points = rng.random((n, d))
    return np.mean(f(points))

dimensions = [1, 2, 3, 5, 10]
budgets = [100, 1000, 10000]
trials = 200
rng = np.random.default_rng(42)

print("Dimension | Budget | Midpoint evals | Midpoint error | MC RMSE")

midpoint_results = {budget: [] for budget in budgets}
mc_results = {budget: [] for budget in budgets}

for d in dimensions:
    exact = d / 3

    for budget in budgets:
        # largest complete midpoint grid that fits within the budget
        m = max(1, int(budget ** (1 / d)))

        while (m + 1)**d <= budget:
            m += 1

        while m**d > budget:
            m -= 1

        midpoint_evals = m**d

        start = time.perf_counter()
        midpoint_estimate = midpoint_rule(d, m)
        midpoint_time = time.perf_counter() - start
        midpoint_error = abs(midpoint_estimate - exact)

        start = time.perf_counter()
        estimates = np.array([
            monte_carlo(d, budget, rng) for _ in range(trials)
        ])
        mc_time = time.perf_counter() - start
        mc_rmse = np.sqrt(np.mean((estimates - exact)**2))

        midpoint_results[budget].append(midpoint_error)
        mc_results[budget].append(mc_rmse)

        print(
            f"{d:9} | {budget:6} | {midpoint_evals:14} | "
            f"{midpoint_error:14.6g} | {mc_rmse:.6g}"
        )

        print(
            f"          Runtime: midpoint={midpoint_time:.4f}s, "
            f"MC ({trials} trials)={mc_time:.4f}s"
        )

for budget in budgets:
    plt.figure(figsize=(8, 5))

    plt.semilogy(
        dimensions,
        midpoint_results[budget],
        "o-",
        label="Midpoint absolute error"
    )
    plt.semilogy(
        dimensions,
        mc_results[budget],
        "s-",
        label="Monte Carlo RMSE"
    )

    plt.xlabel("Dimension")
    plt.ylabel("Error (log scale)")
    plt.title(f"Integration Error vs. Dimension (Budget = {budget:,})")
    plt.xticks(dimensions)
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.show()
