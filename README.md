# PaymentChains

## Overview
This repository accompanies a research project on liquidity frictions in macroeconomics. It provides MATLAB code to solve and simulate a model of payment chains and liquidity shortages. The scripts solve for equilibrium values, perform value-iteration, and analyze continuous-time versions of the model. The code is designed to reproduce results for the associated paper on liquidity frictions and payment-chain dynamics.

## Repository structure
PaymentChains/
├── main_v3_eps.m – baseline equilibrium and transitional dynamics
├── main.m – simple examples of micro objects
├── Val_Iter.m – discrete-time value-iteration driver script
├── Val_Iter_functions.m – helper functions for value iteration
├── Val_Iter_params.m – parameter file for value iteration
├── ContinuousTimeSetting/ – continuous-time model code (HJB, parameters, solutions)
│ ├── main_v1.m, main_v2_production.m – continuous-time simulations
│ ├── CT_solution.m, HJBupdate.m – core continuous-time solvers
│ └── CT_analyticplots.m – analytic plot generation
├── sim_mod.m, sim_planner.m – simulation scripts
├── transitions_v1.m – transitional dynamics
├── production_example.mat, bh_cutoffs_ct.mat – data inputs
└── fig*.pdf/eps – figure outputs


## Requirements
- MATLAB R2020a or later  
- Optimization Toolbox (uses `fsolve`/`fminsearch`)  
- Repository files on MATLAB path

## Running the code

### Basic equilibrium
```matlab
cd PaymentChains
main_v3_eps

## Value iteration
cd PaymentChains
Val_Iter

## Continuous-time model
cd PaymentChains/ContinuousTimeSetting
main_v2_production

Uses CT_parameters.m and CT_solution.m to obtain the continuous-time equilibrium and associated plots.


