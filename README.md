# Integration in High Dimensions

A comparison of two types of numerical integration methods, examining how using these methods to integrate in higher dimensions affects accuracy, convergence, and efficiency.

## Research Question

**How does the dimension of an integration problem affect the relative accuracy and efficiency of midpoint-rule integration and Monte Carlo integration if there is a fixed evaluation budget?**

## Background Information:

Numerical integration becomes much more challenging as the dimension of the integration increases. Deterministic methods, namely the midpoint rule, can achieve very high accuracy in low-dimensional problems, but their computational requirements increase rapidly with dimension. Monte Carlo integration is a very good alternative, as its theoretical convergence rate is independent of dimension, although its accuracy ultimately depends on the variance of the integrand.

This project investigates these differences through mathematical analysis and experiments.
It compares both methods over dimensions 1 through 10, using two multidimensional test functions (sum of squares and product of squares) and several evaluation budgets (100, 1000, 10000).

## Methods

Two integration methods are used:

- **Midpoint Rule:** A _deterministic_ numerical integration method that approximates an integral using function evaluations at the exact centers of subdivisions.
- **Monte Carlo Integration:** A _stochastic_ method that estimates an integral using averages of randomly sampled function values.

The methods are evaluated using two functions over the unit hypercube: a sum of squares and a product of squares.
Both functions have known exact integrals, allowing numerical errors to be calculated directly.

## Experimental Design

The investigation includes:

| Experiment | Configuration |
|---|---|
| Baseline comparison | 1-dimensional integration of \(x^2\) |
| Convergence | Three-dimensional sum of squares |
| Variability | 500 Monte Carlo trials per sample size |
| Multidimensional comparison | Dimensions 1–10, budgets of 100, 1,000, and 10,000 |
| Runtime comparison | Nine timing repetitions per configuration |

Note: the midpoint rule uses the largest possible midpoint grid that does not exceed the specified evaluation budget. Monte Carlo uses the full budget.
Accuracy is also calculated using absolute error, relative error, and root mean square error (RMSE).

## Main Findings

The experiments indicate several important findings:

- In one dimension, the midpoint rule converges much faster than Monte Carlo integration.
- In three dimensions, the convergence rates that were observed very precisely match their theoretical predictions.
- The midpoint rule becomes less viable as the dimension increases, because the number of grid points grows exponentially.
- Monte Carlo maintains a convergence rate proportional to N^(-1/2), independent of dimension. (assuming finite variance).
- The relative performance of both methods always  depends on the integrand. At an evaluation budget of 10,000, the observed crossover occurred between dimensions 4 and 5 for the sum of squares function, and between dimensions 5 and 6 for the product of squares function.

Overall, the results of this project demonstrate that the dimension does not solely determine integration difficulty, as the mathematical properties of the integrand are also crucial.

## Repository Structure

The project is organized as follows:

```text
integration-in-high-dimensions/
├── README.md
├── LICENSE
├── requirements.txt
├── src/
│   ├── monte_carlo.py
│   ├── multidimensional.py
│   ├── multidimensional_product.py
│   ├── convergence_rates.py
│   ├── compare_experiments.py
│   ├── variability_analysis.py
│   └── runtime_comparison.py
├── data/
├── figures/
└── paper/
```

The experiments were implemented in Python using both NumPy and Matplotlib.

To run an experiment:

```bash
git clone https://github.com/jmathew89/integration-in-high-dimensions.git
cd integration-in-high-dimensions
pip install -r requirements.txt
python src/compare_experiments.py
```

## Author

**Jonah Mathew**

[GitHub Repository](https://github.com/jmathew89/integration-in-high-dimensions)
