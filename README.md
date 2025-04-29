# qudit-lc

A complete framework for calculating the Local Clifford (LC) entanglement classes for qutrit graph states up to 7 particles the necessary observables extraction and the corresponding data analysis.

## Overview

The goal of this work was the extraction of the entanglement classes using local complementation for systems with local dimension 3 (qutrits). In principle, there are 3 steps.

- **Orbit extraction:** For this step, we had to calculate the LC orbits and represent them as a graph called Orbit Graph (OG).
- **Define Observables:** At this step, we defined a set of properties emerging from graph theory as well as from quantum theory. For example, the vertex number, the density of the graph, and the Schmidt Measure.
- **Correlation Analysis and Results Presentation:** After the extraction of the datasets, we calculated the Pearson and Kendall's tau coefficients to determine if some of the obtained datasets are correlated. Then some 
Scripts are provided for presenting our results with figures and tables for latex.
- **Benchmarking:** Our code is cross-checked with results already known in the literature for qubits. The results are in agreement with Quantum 4, 305 (2020). 
- **Datasets provided:** All the required datasets are publicly accessible.

## Installation

git clone https://github.com/hrachya-zakaryan/qudit-lc


## Python Version and packages
The user must ensure to use the following versions of the packages.
Python                    3.12.6
matplotlib                3.7.3
networkx                  3.2.1
numpy                     1.25.2
pandas                    2.2.3
scipy                     1.11.2
sinter                    1.14.0
sympy                     1.12
wolframclient             1.4.0

qudit-lc/
│
├── og_data/              <- Data produced for the statistical analysis in txt files.
├── c/               <- File with code in C to speed up the code.
├── orbitsd{d}_n{n}_separated/ <- These files include all orbits for d the local dimension and n the number of particles.  
    ├── orbit_i <- The index i indicates the different orbtis.         
├── aspl.ipynb <- Notebook for the extraction of the average shortest path and the heuristic method adopted.
├── data_analysis_functions.py <- A module with every function required for the data extraction and analysis.
├── main_data_notebook.ipynb <- The notebook with the statistical analysis after the data extraction.
├── min_edge_color_states_wolfram_engine.py <- Module to compute the minimum edge chromatic number.
├── mspl.ipynb <- Notebook for the extraction of the maximun shortest path and the heuristic method adopted.
├── real_join_data.py <- Module for the extraction of the required datasets.
└── sm_final_data.py <- Module for preprocessing the Schmidt Measure data before the data analysis.

