# Analysing-Plasma-Channels

## Abel-Inversion.py

This code performs an Abel inversion on a vertical phase strip for a plasma channel using the PyAbel package. The input data (i.e. the phase profile) need to be a .txt file. Smoothing can be applied to the phase and a simple Monte Carlo uncertainty analysis can be performed.

## Plasma-columnn.ipynb

This code contains 2 labelled parts. The first part plots a 3D surface plot of the normalised intensity of a Gaussian Beam.
The second part plots the 2D shape of a Helium plasma column generated through barrier suppression. When plotting, a slider widget is available to easily change the inensity of the gaussian beam to see the changes in the shape of the plasma column. The data is calculated for ony 1 Rayleigh range but can be changed if need be.
