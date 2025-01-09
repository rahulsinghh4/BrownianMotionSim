# Brownian Motion Simulation

This MATLAB-based simulation models the Brownian motion of particles suspended in a fluid, designed to complement experimental observations of polystyrene beads in glycol solution. The simulation provides theoretical predictions that can be compared with experimental data gathered using computer vision tracking.

## Overview

The simulation implements both 1D and 2D Brownian motion models, incorporating various physical parameters and statistical analyses. It's particularly useful for validating experimental setups and understanding the underlying physics of particle diffusion.

## Physical Parameters

- Particle diameter (d): 1.0e-6 meters
- Fluid viscosity (η): 1.0e-3 Pascal-seconds (water at 293K)
- Temperature (T): 293 Kelvin
- Boltzmann constant (kB): 1.38e-23 J/K

The theoretical diffusion coefficient (D) is calculated using the Stokes-Einstein relation:
```matlab
D = kB * T / (3 * π * η * d)
```

## Key Features

### 1. Basic Simulations
- One-dimensional Brownian motion with displacement plotting
- Two-dimensional particle tracking
- Squared displacement calculations and analysis

### 2. Advanced Features
- Multiple particle simulations (default: 10 particles)
- Bulk flow effects simulation
- Statistical analysis including:
  - Mean squared displacement
  - Ensemble averaging
  - Error estimation
  - Auto-correlation analysis
  - Cross-correlation between particles

### 3. Statistical Tools
- Sampling uncertainty visualization
- Chi-squared distribution analysis
- Population mean studies
- Correlation analysis for particle trajectories

## Implementation Details

### Single Particle Simulation
The simulation uses random normal distributions to generate particle displacements:
```matlab
dx = k * randn(N,1)
dy = k * randn(N,1)
```
where `k = sqrt(D * dimensions * τ)`, and τ is the time interval.

### Multiple Particle Analysis
The code implements a cell array structure to track multiple particles:
```matlab
particle{i}.dx = k * randn(1,N)
particle{i}.x = cumsum(particle{i}.dx)
```

### Error Analysis
Several error metrics are calculated:
- Standard error in diffusion coefficient
- Actual error compared to theoretical value
- Ensemble averaging across multiple particles

## Experimental Validation

This simulation is designed to be compared with experimental data from:
- Tracking system: C#/.NET-based computer vision program
- Particles: Polystyrene beads
- Medium: Glycol solution
- Measurement: Real-time particle position tracking

## Advanced Analysis Tools

### Correlation Analysis
- Auto-correlation functions for single particle trajectories
- Cross-correlation between different particles
- Correlated trajectory generation with adjustable correlation coefficient

### Statistical Distribution
- Histogram analysis of particle displacements
- Population mean studies
- Sampling uncertainty visualization

## Output Visualizations

The simulation provides several visualization options:
1. Particle trajectories (1D and 2D)
2. Displacement squared vs. time
3. Histogram distributions
4. Correlation function plots
5. Error analysis plots

## Notes on Experimental Comparison

When comparing with experimental data:
1. Ensure matching time intervals between simulation and experiment
2. Account for camera frame rate in experimental setup
3. Consider bulk flow effects in the experimental medium
4. Compare statistical distributions rather than individual trajectories
5. Use error bars for meaningful comparison

## Limitations

- Assumes perfect spherical particles
- Neglects particle-particle interactions
- Simplified fluid dynamics
- Ideal temperature distribution
- No boundary effects considered

## Future Improvements

Possible enhancements to the simulation:
1. Implementation of particle-particle interactions
2. Addition of boundary conditions
3. Temperature gradient effects
4. Non-spherical particle dynamics
5. Multiple particle size distributions

## Dependencies

- MATLAB (core functionality)
- Statistics and Machine Learning Toolbox (for advanced statistical analysis) in MATLAB
