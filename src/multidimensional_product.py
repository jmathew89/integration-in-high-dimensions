import numpy as np
import matplotlib.pyplot as plt
import time

def g(x):
    return np.prod(x**2, axis=1)

def midpoint_rule(d, m):
    coords = (np.arange(m) + 0.5) / m
    grids = np.meshgrid(*([coords] * d), indexing="ij")
    points = np.stack([grid.ravel() for grid in grids], axis=1)
    return np.mean(g(points))

def monte_carlo(d, n, rng):
    points = rng.random((n, d))
    return np.mean(g(points))

dimensions = [1, 2, 3, 5, 10]
budgets = [100, 1000, 10000]
trials = 200
rng = np.random.default_rng(42)

midpoint_abs = {b: [] for b in budgets}
mc_abs = {b: [] for b in budgets}
midpoint_rel = {b: [] for b in budgets}
mc_rel = {b: [] for b in budgets}

print(
    "Dimension | Budget | Midpoint evals | "
    "Midpoint abs. error | MC RMSE | Midpoint rel. error | MC relative RMSE"
)

for d in dimensions:
    exact = (1 / 3)**d

    for budget in budgets:
        # largest complete grid that fits within the budget
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
        midpoint_relative = midpoint_error / exact

        start = time.perf_counter()
        estimates = np.array([
            monte_carlo(d, budget, rng)
            for _ in range(trials)
        ])
        mc_time = time.perf_counter() - start

        errors = estimates - exact
        mc_rmse_value = np.sqrt(np.mean(errors**2))
        mc_relative_value = mc_rmse_value / exact

        midpoint_abs[budget].append(midpoint_error)
        mc_abs[budget].append(mc_rmse_value)
        midpoint_rel[budget].append(midpoint_relative)
        mc_rel[budget].append(mc_relative_value)

        print(
            f"{d:9} | {budget:6} | {midpoint_evals:14} | "
            f"{midpoint_error:.6g} | {mc_rmse_value:.6g} | "
            f"{midpoint_relative:.6g} | {mc_relative_value:.6g}"
        )

        print(
            f"          Runtime: midpoint={midpoint_time:.4f}s, "
            f"MC ({trials} trials)={mc_time:.4f}s"
        )

# plot relative error.
for budget in budgets:
    plt.figure(figsize=(8, 5))

    plt.semilogy(
        dimensions,
        midpoint_rel[budget],
        "o-",
        label="Midpoint relative error"
    )
    plt.semilogy(
        dimensions,
        mc_rel[budget],
        "s-",
        label="Monte Carlo relative RMSE"
    )

    plt.xlabel("Dimension")
    plt.ylabel("Relative error (log scale)")
    plt.title(f"Product Function: Relative Error (Budget = {budget:,})")
    plt.xticks(dimensions)
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.show()
    plt.close()
