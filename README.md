# MATLAB Script for Coherence and PSD Analysis

## Overview
This MATLAB script processes mechanomyography (MMG) data from multiple participants to analyze coherence and power spectral density (PSD) metrics. The script performs data preprocessing, segmentation, and analysis of muscle activation patterns for different gestures, including statistical and coherence-based evaluations.

The entire code in this repository is based on the PhD work of Carlos Sebastian Mancero Castillo.

## Features
- **Data Preprocessing:**
  - Loads MMG data from different participants.
  - Filters data using bandpass filtering (1-150 Hz).
  - Segments data based on predefined gestures and force levels.

- **Signal Processing:**
  - Computes Hilbert transform for envelope extraction.
  - Applies Welch’s method for spectral estimation.
  - Performs FFT-based power spectral density analysis.

- **Coherence Analysis:**
  - Computes magnitude-squared coherence (MSC) between muscle pairs.
  - Uses confidence intervals to validate coherence results.
  - Implements non-negative matrix factorization (NNMF) to extract signal components.

- **Statistical and Graphical Analysis:**
  - Generates coherence heatmaps and connectivity matrices.
  - Produces box plots for statistical comparisons across conditions.
  - Uses adjacency matrices for network connectivity analysis.

## Data Processing Steps
1. **Load and Preprocess Data:**
   - Loads MMG signals from CSV files.
   - Filters and segments data into trials.
   
2. **Compute Coherence and PSD:**
   - Calculates coherence for muscle pairs.
   - Applies Welch’s method for spectral density estimation.
   
3. **Statistical Analysis:**
   - Computes confidence intervals.
   - Performs network connectivity analysis with betweenness centrality and clustering coefficients.
   
4. **Plotting and Visualization:**
   - Generates coherence plots.
   - Produces heatmaps for adjacency matrices.
   - Compares statistical differences using box plots.

## Dependencies
- MATLAB Signal Processing Toolbox
- Butterworth filter functions
- Tukey windowing for signal tapering
- Non-negative matrix factorization (NNMF)

## File Structure
- **Participants/**: Contains individual participant data files.
- **Data Library/**: Stores raw MMG datasets.
- **Results/**: Stores computed coherence and PSD results.
- **Figures/**: Saves generated plots and visualizations.

## How to Run the Script
1. Ensure that all required data files are present in the appropriate directories.
2. Modify participant selection (`who` variable) as needed.
3. Run the script in MATLAB.
4. Processed data and figures will be saved automatically.

## Output Files
- `*-coh.mat`: Coherence results per participant.
- `*-psd.mat`: PSD results per participant.
- `*-params.mat`: Parameters used in analysis.
- `*-nf.mat`: Frequency vector for coherence computation.

## Acknowledgments
This script was developed for the analysis of MMG coherence and PSD for connectivity and network analysis.

