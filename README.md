# Bayesian Parameter Estimation of Network SIS Model

>*Note: This README is AI generated!

This project implements Bayesian parameter estimation for a Network Susceptible-Infected-Susceptible (SIS) epidemic model using Markov Chain Monte Carlo (MCMC) methods and an Extended Kalman Filter (EKF). The aim is to estimate critical parameters such as transmission and recovery rates within a network-based SIS model framework.

## Project Structure

- **ekf/**: Contains MATLAB implementation for the Extended Kalman Filter method for online parameter and state estimation.
- **mcmc/**: Holds MATLAB scripts for Bayesian parameter estimation using Markov Chain Monte Carlo methods.
- **mcmc_seir/**: Related extensions or exploratory implementations involving SEIR-type epidemic models.
- **papers/**: Collection of relevant research papers and reference materials.
- **report/**: LaTeX and PDF documents detailing comprehensive results, simulations, and comparative analyses.

## Key Components

### Network SIS Model
- Investigates disease spread dynamics on interconnected networks.
- Uses differential equations capturing infection and recovery dynamics.

### Bayesian Estimation with Markov Chain Monte Carlo (MCMC)
- Implements Random-Walk Metropolis-Hastings algorithm.
- Provides posterior distributions to quantify parameter uncertainty.

### Extended Kalman Filter (EKF)
- Performs sequential parameter and state estimation from noisy observations.
- Includes calculation of system and measurement Jacobians.

