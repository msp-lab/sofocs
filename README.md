# Sparse Optimization for Fractional-Order Chaotic System Implementation (SOFOCS)

## Background
This is Matlab implementation for paper "Analysis and Implementation of Fractional-Order Chaotic System with Standard Components". The paper has been submitted to the Journal of Advanced Research for review.

The new approach is very helpful for implementing arbitrary fractional chaotic systems with commercially available components. Moreover, the system accuracy and complexity can be analyzed by simulation.


## Function Description

- sofocs.m: Finding an Approximation Optimal Solution for the Franctance in Chaotic Circuit System

- mag_curve.m:  Transfer Function H(s) Curve for the Fractional-order Fractance

- objective.m: Objective Function

- uncertainty.m: Monte Carlo based Unvertainty Estimation Algorithm

- tf_zpk.m: Transfer Function H(s) of the "pole/zero" Method

- value2code.m: Convert Component Values to Parameter Matrix X

- rtnorm.m: Pseudorandom numbers from a truncated Gaussian distribution

- optim_plot.m: Plot the result of optimal values


## Quick Start

### Setting Parameters

1. Change the parameters in main.m

	q = 0.75;    % fractional order
	
	delta = 1;   % dB, The maximum discrepancy
	
	structure = 'Chain'; % Fractance structure {'Chain', 'Tree', 'Ladder'}
	
	N = 4;       % system order of franctance

2. Change the threhold of uncertainty

	MaxUC = 0.2;    

### Calculating

	Run main.m for a quick example.

## Contact Author

- Kunpeng Wang, e-mail: wkphnzk@163.com

