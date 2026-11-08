import numpy as np
import matplotlib.pyplot as plt

# experiment settings
d = 3
sample_sizes = [100, 300, 1000, 3000, 10000]
trials = 500

rng = np.random.default_rng(42)


def sum_squares(x):
    return np.sum(x**2, axis=1)


exact = d / 3


def monte_carlo(n):
    points = rng.random((n, d))
    return np.mean(sum_squares(points))


# Store results
observed_std = []
theoretical_std = []

print("Sample size | Observed SD | Theoretical SD")
print("-" * 48)

for n in sample_sizes:
    estimates = np.array([
        monte_carlo(n)
        for _ in range(trials)
    ])

    # st.dev across independent trials
    obs = np.std(estimates, ddof=1)

    # Theoretical st. dev
    theory = np.sqrt(4 * d / (45 * n))

    observed_std.append(obs)
    theoretical_std.append(theory)

    print(f"{n:11} | {obs:.6f} | {theory:.6f}")

    # estimate convergence rate
observed_slope = np.polyfit(
    np.log(sample_sizes), np.log(observed_std), 1
)[0]

theoretical_slope = np.polyfit(
    np.log(sample_sizes), np.log(theoretical_std), 1
)[0]

print(f"\nObserved log-log slope:    {observed_slope:.4f}")
print(f"Theoretical log-log slope: {theoretical_slope:.4f}")
print("Expected theoretical slope: -0.5000")

plt.figure(figsize=(8, 5))
plt.loglog(sample_sizes, observed_std, "o-", label="Observed SD")
plt.loglog(sample_sizes, theoretical_std, "s--", label="Theoretical SD")

plt.xlabel("Sample size N")
plt.ylabel("Standard deviation")
plt.title("Monte Carlo Variability vs. Sample Size (d = 3)")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
