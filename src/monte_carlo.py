import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**2

def midpoint_rule(n):
    dx = 1 / n
    x = (np.arange(n) + 0.5) * dx
    return dx * np.sum(f(x))

def monte_carlo(n, rng):
    x = rng.uniform(0, 1, size=n)
    return np.mean(f(x))

exact = 1 / 3
sample_sizes = np.array([10, 30, 100, 300, 1000, 3000, 10000])
trials = 200
rng = np.random.default_rng(42)

midpoint_errors = []
mc_rmse = []

for n in sample_sizes:
    midpoint_estimate = midpoint_rule(n)
    midpoint_errors.append(abs(midpoint_estimate - exact))

    estimates = np.array([
        monte_carlo(n, rng) for _ in range(trials)
    ])

    errors = estimates - exact
    mc_rmse.append(np.sqrt(np.mean(errors**2)))

plt.figure(figsize=(8, 5))
plt.loglog(sample_sizes, midpoint_errors, "o-", label="Midpoint error")
plt.loglog(sample_sizes, mc_rmse, "o-", label="Monte Carlo RMSE")

plt.xlabel("Function evaluations (N)")
plt.ylabel("Error")
plt.title("Integration Error vs. Computational Cost")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.show()

print("N       Midpoint error       Monte Carlo RMSE")
for n, m_error, mc_error in zip(
    sample_sizes, midpoint_errors, mc_rmse
):
    print(f"{n:<7} {m_error:<20.8g} {mc_error:.8g}")
