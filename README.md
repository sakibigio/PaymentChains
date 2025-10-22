# Payment Chains Model

**Authors:** Saki Bigio (UCLA), Esteban Méndez (Central Bank of Costa Rica), Diana Van Patten (Yale University)

## Overview

This repository contains MATLAB code to replicate the results from **"A Theory of Payments-Chain Crises"**. The paper introduces an endogenous network of payment chains into a business cycle model, analyzing how delayed payments create production delays and aggregate productivity losses.

The code solves for equilibrium allocations using value function iteration in both endowment and production economies, computes transitional dynamics, and characterizes the optimal Ramsey planner solution.

## Repository Structure

```
PaymentChains/
├── run_entrepreneur.m              # Main script: Endowment economy (value iteration)
├── run_entrepreneur_production.m   # Main script: Production economy (value iteration)
├── run_transitions.m               # Main script: Transitional dynamics
├── run_planner.m                   # Main script: Ramsey planner solution
│
├── functions/                      # Helper functions and parameters
│   ├── Val_Iter_params.m          # Parameters for value iteration
│   ├── Val_Iter_functions.m       # Core functions for value iteration
│   ├── Val_iter_sim.m             # Simulation functions
│   ├── Trans_params.m             # Parameters for transitions
│   ├── Trans_function.m           # Transition helper functions
│   └── [other utility functions]
│
├── plotting/                       # Plotting functions
│   ├── Val_Iter_plots.m           # Value iteration plots
│   ├── Val_Iter_simplots.m        # Simulation plots
│   ├── Plots_AggregateEuler.m     # Aggregate Euler equation plots
│   └── [other plotting utilities]
│
├── tests/                          # Testing and validation
│   └── Val_Iter_tests.m
│
├── data/                           # Input data
│   └── bh_cutoffs_ct.mat          # Cutoffs from continuous-time analysis
│
├── archived/                       # Archived code (not for publication)
│   ├── ContinuousTimeSetting/     # Continuous-time model (testing only)
│   ├── Old Codes/
│   └── Old Figures/
│
├── Matlab Figures/                 # Output directory for generated figures
└── production_example.mat          # Example production data
```

## Requirements

- **MATLAB** R2020a or later
- **Optimization Toolbox** (uses `fsolve`, `fminunc`, `fmincon`)
- All repository folders added to MATLAB path (handled automatically by scripts)

## Quick Start

1. Clone or download this repository
2. Open MATLAB and navigate to the repository root directory
3. Run any of the main scripts:

### Endowment Economy (Value Iteration)
```matlab
run_entrepreneur
```
Solves the household problem with endowments using value function iteration. Computes policy functions and value functions for debt choices.

**Key outputs:**
- Value functions and policy functions
- Equilibrium thresholds (B*, B̃)
- Simulated household trajectories

### Production Economy (Value Iteration)
```matlab
run_entrepreneur_production
```
Extends the model to include capital investment decisions. Solves the two-dimensional value function over debt and capital.

**Key outputs:**
- 2D value function V(B,K)
- Joint policy functions for debt and capital
- Comparative statics

### Transitional Dynamics
```matlab
run_transitions
```
Computes the economy's transitional dynamics following shocks to the borrowing limit or other parameters.

**Key outputs:**
- Impulse response functions
- Phase diagrams
- Transition paths for key aggregates (Y, C, B, μ)

### Ramsey Planner Solution
```matlab
run_planner
```
Solves for the constrained-efficient allocation chosen by a Ramsey planner who internalizes payment-chain externalities.

**Key outputs:**
- Optimal debt policy B*(B̃)
- Comparison of decentralized vs. planner solutions
- Welfare analysis

## Model Overview

The model features:
- **Payment chains**: Firms can delay payments until they receive payment from downstream customers
- **Endogenous TFP**: The fraction of chained payments (μ) determines aggregate productivity through A(μ)
- **Liquidity constraints**: Limited ability to borrow creates trade-offs in payment timing
- **Externalities**: Private agents don't internalize the costs their payment delays impose on others

**Key equations:**
- TFP function: A(μ) = δ(1-μ)/(1-δμ)
- Chained payment fraction: μ = 1 - (1-β)B - S^w
- Natural borrowing limit: B̄ = 1/(1-β)

## Key Parameters

Configurable in `functions/Val_Iter_params.m` and `functions/Trans_params.m`:

- `beta`: Household discount factor (default: 0.95)
- `delta`: Production discount factor for chained payments (default: 0.9)
- `theta`: Pareto weight on workers (default: 0.75)
- `alpha`: Capital share in production (default: 0.5)
- `q`: Price of chained vs. spot payments (default: 1.25)

## Output

Generated figures are saved to `Matlab Figures/` and include:
- Policy functions (debt, capital, consumption)
- Value functions
- Transition dynamics
- Comparative statics
- Phase diagrams

Figures can be automatically synced to Overleaf using `utils/Plot_Migration.m` (edit paths as needed).

## Citation

If you use this code in your research, please cite:

```
Bigio, S., Méndez, E., & Van Patten, D. (2024).
A Theory of Payments-Chain Crises.
Working Paper.
```

## Contact

- **Saki Bigio**: sbigio@econ.ucla.edu (UCLA)
- **Esteban Méndez**: mendezce@bccr.fi.cr (Central Bank of Costa Rica)
- **Diana Van Patten**: diana.vanpatten@yale.edu (Yale University)

## Notes

- The `archived/ContinuousTimeSetting/` folder contains continuous-time versions of the model used for testing and robustness checks. These are not required for replicating the main results.
- All main scripts automatically add necessary folders to the MATLAB path
- Computation times vary: endowment economy (~1-2 minutes), production economy (~5-10 minutes depending on grid resolution)

## License

This code is provided for academic research purposes. Please contact the authors for questions about usage rights.

---

**Version:** October 2024  
**Repository:** https://github.com/sakibigio/PaymentChains
