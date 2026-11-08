import numpy as np
import matplotlib.pyplot as plt

# dimensions and evaluation budgets for the experiments
dimensions = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
budgets = [100, 1000, 10000]
trials = 300

rng = np.random.default_rng(42)


# sum of squares
# Exact integral over [0, 1]^d is d / 3
def sum_squares(x):
    return np.sum(x**2, axis=1)


def sum_squares_exact(d):
    return d / 3


# product of squares
# Exact integral over [0, 1]^d is (1 / 3)^d
def product_squares(x):
    return np.prod(x**2, axis=1)


def product_squares_exact(d):
    return (1 / 3) ** d

def midpoint_rule(function, d, m):
    """Approximate an integral using m midpoint subdivisions per dimension."""
    coords = (np.arange(m) + 0.5) / m

    grids = np.meshgrid(*([coords] * d), indexing="ij")
    points = np.stack([grid.ravel() for grid in grids], axis=1)

    return np.mean(function(points))


def monte_carlo(function, d, n, rng):
    """Approximate an integral using n uniform random samples."""
    points = rng.random((n, d))
    return np.mean(function(points))

# store results from each experiment
results = []

test_functions = [
    ("Sum of squares", sum_squares, sum_squares_exact),
    ("Product of squares", product_squares, product_squares_exact),
]

for function_name, function, exact_function in test_functions:
    for d in dimensions:
        exact = exact_function(d)

        for budget in budgets:
            # chooses the largest complete midpoint grid within the budget.
            m = int(budget ** (1 / d))

            # float rounding near integer roots.
            while (m + 1) ** d <= budget:
                m += 1
            while m**d > budget:
                m -= 1

            midpoint_n = m**d
            midpoint_estimate = midpoint_rule(function, d, m)
            midpoint_abs_error = abs(midpoint_estimate - exact)
            midpoint_rel_error = midpoint_abs_error / abs(exact)

            # Run monte carlo independently for each trial.
            mc_estimates = np.array([
                monte_carlo(function, d, budget, rng)
                for _ in range(trials)
            ])

            mc_rmse = np.sqrt(np.mean((mc_estimates - exact) ** 2))
            mc_rel_rmse = mc_rmse / abs(exact)

            # variability across independent Monte Carlo estimates
            mc_std = np.std(mc_estimates, ddof=1)
            mc_mean = np.mean(mc_estimates)

            mc_bias = mc_mean - exact

            mc_ci_low = np.percentile(mc_estimates, 2.5)
            mc_ci_high = np.percentile(mc_estimates, 97.5)

            results.append({
                "function": function_name,
                "dimension": d,
                "budget": budget,
                "midpoint_evaluations": midpoint_n,
                "midpoint_abs_error": midpoint_abs_error,
                "midpoint_rel_error": midpoint_rel_error,
                "mc_evaluations_per_trial": budget,
                "mc_rmse": mc_rmse,
                "mc_rel_rmse": mc_rel_rmse,
                "mc_std": mc_std,
                "mc_bias": mc_bias,
                "mc_trial_2.5_percentile": mc_ci_low,
                "mc_trial_97.5_percentile": mc_ci_high,
            })

print("Experiments completed:", len(results))

# Print a summary
print("\nSUMMARY OF RESULTS")
print("=" * 105)

for result in results:
    print(
        f"{result['function']:18} | "
        f"d={result['dimension']:2} | "
        f"budget={result['budget']:5} | "
        f"midpoint n={result['midpoint_evaluations']:5} | "
        f"mid abs err={result['midpoint_abs_error']:.6g} | "
        f"mid rel err={result['midpoint_rel_error']:.6g} | "
        f"MC RMSE={result['mc_rmse']:.6g} | "
        f"MC rel RMSE={result['mc_rel_rmse']:.6g}"
        f" | MC std={result['mc_std']:.6g}"
        f" | MC bias={result['mc_bias']:.6g}"
        f" | MC 95% range=["
        f"{result['mc_trial_2.5_percentile']:.6g}, "
        f"{result['mc_trial_97.5_percentile']:.6g}]"
    )

    # plot both functions
plot_budget = 10000

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

plot_functions = [
    ("Sum of squares", 0),
    ("Product of squares", 1),
]

for row, (function_name, _) in enumerate(plot_functions):
    selected = [
        r for r in results
        if r["function"] == function_name
        and r["budget"] == plot_budget
    ]

    dims = [r["dimension"] for r in selected]

    midpoint_abs = [r["midpoint_abs_error"] for r in selected]
    mc_abs = [r["mc_rmse"] for r in selected]

    midpoint_rel = [r["midpoint_rel_error"] for r in selected]
    mc_rel = [r["mc_rel_rmse"] for r in selected]

    # Absolute error
    axes[row, 0].plot(
        dims, midpoint_abs, "o-", label="Midpoint"
    )
    axes[row, 0].plot(
        dims, mc_abs, "s-", label="Monte Carlo RMSE"
    )
    axes[row, 0].set_title(f"{function_name}: Absolute Error")
    axes[row, 0].set_xlabel("Dimension")
    axes[row, 0].set_ylabel("Error")
    axes[row, 0].set_yscale("log")
    axes[row, 0].legend()
    axes[row, 0].grid(True, which="both", alpha=0.3)

    # Relative error
    axes[row, 1].plot(
        dims, midpoint_rel, "o-", label="Midpoint"
    )
    axes[row, 1].plot(
        dims, mc_rel, "s-", label="Monte Carlo relative RMSE"
    )
    axes[row, 1].set_title(f"{function_name}: Relative Error")
    axes[row, 1].set_xlabel("Dimension")
    axes[row, 1].set_ylabel("Relative error")
    axes[row, 1].set_yscale("log")
    axes[row, 1].legend()
    axes[row, 1].grid(True, which="both", alpha=0.3)

fig.suptitle(
    f"Integration Error vs. Dimension (Budget = {plot_budget})",
    fontsize=14
)
fig.tight_layout()
plt.show()

# find crossover ranges with the largest evaluation budget
print("\nCROSSOVER ANALYSIS")
print("=" * 65)

for function_name, _, _ in test_functions:
    selected = sorted(
        [
            r for r in results
            if r["function"] == function_name
            and r["budget"] == plot_budget
        ],
        key=lambda r: r["dimension"]
    )

    print(f"\n{function_name}:")

    previous_winner = None

    for r in selected:
        midpoint_error = r["midpoint_abs_error"]
        mc_error = r["mc_rmse"]

        winner = (
            "Monte Carlo"
            if mc_error < midpoint_error
            else "Midpoint"
        )

        print(
            f"d={r['dimension']:2}: "
            f"midpoint={midpoint_error:.6g}, "
            f"MC RMSE={mc_error:.6g} -> {winner}"
        )

        if previous_winner is not None and winner != previous_winner:
            print(
                f"  ** Winner changes between dimensions "
                f"{previous_dimension} and {r['dimension']} **"
            )

        previous_winner = winner
        previous_dimension = r["dimension"]
