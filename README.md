# ATP-synthase

This repository contains all of the code which generates results for the MITACS GRI 2023 project "Design principles of molecular machines: Optimal gearing ratio for coupled rotary motors" by J.C. Perez-Ramirez, M.P. Leighton and David A. Sivak, based on the code which generates results for the manuscript titled "Nonequilibrium Response of Stochastic Strongly-Coupled Rotary Motors." by E. Lathouwers , J.N.E. Lucero, and D.A. Sivak.

# Directories reference:

src: directory containing the scripts used; in EcoupleValue the ones that implement the integration of the Smoluchowski equation (using FTCS) and the associated analysis discussed in the manuscript. 

graphs: directory containing the graphs produced during the project, one for the work outflow and efficiency of ATP synthase given different Ecouple values and the tight coupling limit curve.

# Description

A separate EcoupleValue directory is expected for each Ecouple value tested.

main.py calculates the near steady state density function for the system using the fpe and utilities modules (once compiled). 

plot_ATP_ratio.py processes the density function file to calculate the flux, which is used to get the power outflow and the efficiency of the system for each of the n0 values in the array, all as described in the manuscript. Finally, it plots the power and efficiency curves for the given parameters.

plotEcouples.py plots the curves for all of the different Ecouple values in the same graph: it needs the data produced by plot_ATP_ratio.py.

tightCoupling.py calculates and plots the power and efficiency curves for the tight coupling limit of the system given the different n0 values.

Each of the bash run files submits the referenced script to the cluster.

