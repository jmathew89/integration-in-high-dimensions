import time
import numpy as np
import matplotlib.pyplot as plt

dimensions = list(range(1, 11))
budgets = [100, 1000, 10000]
repeats = 9

rng = np.random.default_rng(2026)


def sum_squares(x):
    return np.sum(x**2, axis=1)


def midpoint_rule(d, m):
    coords = (np.arange(m) + 0.5) / m
    grids = np.meshgrid(*([coords] * d), indexing="ij")
    points = np.stack([grid.ravel() for grid in grids], axis=1)
    return np.mean(sum_squares(points))


def monte_carlo(d, n):
    points = rng.random((n, d))
    return np.mean(sum_squares(points))



def median_runtime(function):
    function()  # Warm-up

    batch_size = 50
    times = []

    for _ in range(repeats):
        start = time.perf_counter()

        for _ in range(batch_size):
            function()

        elapsed = time.perf_counter() - start
        times.append(elapsed / batch_size)

    return np.median(times)



results = []

print(
    f"{'d':>2} {'Budget':>7} {'Midpoint evals':>15} "
    f"{'Midpoint (ms)':>15} {'MC (ms)':>12}"
)
print("-" * 58)

for d in dimensions:
    for budget in budgets:
        # Largest possible grid that stays within budget
        m = 1
        while (m + 1) ** d <= budget:
            m += 1

        midpoint_evals = m**d

        midpoint_time = median_runtime(
            lambda: midpoint_rule(d, m)
        )

        mc_time = median_runtime(
            lambda: monte_carlo(d, budget)
        )

        results.append({
            "dimension": d,
            "budget": budget,
            "midpoint_evals": midpoint_evals,
            "midpoint_ms": midpoint_time * 1000,
            "mc_ms": mc_time * 1000,
        })

        print(
            f"{d:2d} {budget:7d} {midpoint_evals:15d} "
            f"{midpoint_time * 1000:15.3f} "
            f"{mc_time * 1000:12.3f}"
        )


# plot runtime by dimension, one panel for each evaluation budget (100, 1,000, 10,000). 
# log scale for yaxis
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for ax, budget in zip(axes, budgets):
    subset = [r for r in results if r["budget"] == budget]

    ax.plot(
        [r["dimension"] for r in subset],
        [r["midpoint_ms"] for r in subset],
        "o-", label="Midpoint"
    )
    ax.plot(
        [r["dimension"] for r in subset],
        [r["mc_ms"] for r in subset],
        "s-", label="Monte Carlo"
    )

    ax.set_title(f"Budget = {budget:,}")
    ax.set_xlabel("Dimension")
    ax.set_ylabel("Median runtime (ms)")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend()

plt.tight_layout()
plt.show()
