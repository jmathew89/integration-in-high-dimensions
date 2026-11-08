import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sum(x**2, axis=1)

def midpoint_rule(d, m):
    coords = (np.arange(m) + 0.5) / m
    grids = np.meshgrid(*([coords] * d), indexing="ij")
    points = np.stack([grid.ravel() for grid in grids], axis=1)
    return np.mean(f(points))

def monte_carlo(d, n, rng):
    points = rng.random((n, d))
    return np.mean(f(points))

d = 3
exact = d / 3
budgets = np.array([27, 64, 125, 216, 343, 512, 729, 1000, 1728, 4096, 8000])
trials = 300
rng = np.random.default_rng(42)

midpoint_n = []
midpoint_errors = []
mc_rmse = []

for budget in budgets:
    # use a complete grid with m ppd.
    m = int(round(budget ** (1 / d)))
    actual_n = m**d

    midpoint_estimate = midpoint_rule(d, m)
    midpoint_error = abs(midpoint_estimate - exact)

    estimates = np.array([
        monte_carlo(d, actual_n, rng)
        for _ in range(trials)
    ])
    rmse = np.sqrt(np.mean((estimates - exact)**2))

    midpoint_n.append(actual_n)
    midpoint_errors.append(midpoint_error)
    mc_rmse.append(rmse)

log_n = np.log(np.array(midpoint_n))
mid_slope = np.polyfit(log_n, np.log(midpoint_errors), 1)[0]
mc_slope = np.polyfit(log_n, np.log(mc_rmse), 1)[0]

print(f"\nObserved midpoint slope: {mid_slope:.4f}")
print(f"Observed Monte Carlo slope: {mc_slope:.4f}")
print("Theoretical midpoint slope: -0.6667")
print("Theoretical Monte Carlo slope: -0.5000")

print("Evaluations | Midpoint error | Monte Carlo RMSE")
for n, me, mc in zip(midpoint_n, midpoint_errors, mc_rmse):
    print(f"{n:11} | {me:.6g} | {mc:.6g}")

plt.figure(figsize=(8, 5))
plt.loglog(midpoint_n, midpoint_errors, "o-", label="Midpoint error")
plt.loglog(midpoint_n, mc_rmse, "s-", label="Monte Carlo RMSE")

n_ref = np.array([min(midpoint_n), max(midpoint_n)])
mid_ref = (1 / (4 * midpoint_n[0]**(2/d))) * (
    n_ref / midpoint_n[0]
)**(-2/d)

mc_ref = np.sqrt(4 * d / (45 * midpoint_n[0])) * (
    n_ref / midpoint_n[0]
)**(-0.5)

plt.loglog(n_ref, mid_ref, "--", label="Reference slope: N^(-2/3)")
plt.loglog(n_ref, mc_ref, "--", label="Reference slope: N^(-1/2)")

plt.xlabel("Actual function evaluations")
plt.ylabel("Error (log scale)")
plt.title("Convergence Rates in 3 Dimensions")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.show()
