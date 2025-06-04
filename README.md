# qudit-lc

A complete framework for calculating the Local Clifford (LC) entanglement classes for qudit graph states, the necessary observables extraction, and the corresponding data analysis.

## Overview

The goal of this work was the extraction of the entanglement classes using local complementation. In principle, there are 3 steps.

- **Orbit extraction:** For this step, we had to calculate the LC orbits and represent them as a graph called Orbit Graph (OG).
- **Define Observables:** At this step, we defined a set of properties emerging from graph theory as well as from quantum theory. For example, the vertex number, the density of the graph, and the Schmidt Measure.
- **Correlation Analysis and Results Presentation:** After the extraction of the datasets, we calculated the Pearson and Kendall's tau coefficients to determine if some of the obtained datasets are correlated. Then some 
Scripts are provided for presenting our results with figures and tables for latex.
- **Benchmarking:** Our code is cross-checked with results already known in the literature for qubits. The results are in agreement with Quantum 4, 305 (2020). 
- **Datasets provided:** All the required datasets are publicly accessible.

## Installation

```git clone https://github.com/hrachya-zakaryan/qudit-lc```


## Python Version and Packages
The user must ensure to use the following versions of the packages.<br>
Python                    3.11.5 <br>
matplotlib                3.7.3 <br>
networkx                  3.2.1 <br>
numpy                     1.25.2 <br>
pandas                    2.2.3 <br>
scipy                     1.11.2 <br>
sinter                    1.14.0 <br>
sympy                     1.12 <br>
wolframclient             1.4.0 <br>

## File Hierarchy
qudit-lc/ <br>
&emsp;&emsp;│ <br>
&emsp;&emsp;├── og_data_final/   <- Data produced for the statistical analysis in txt files. <br>
&emsp;&emsp;├── c/               <- Folder with code in C for graph generation and complementation. <br>
&emsp;&emsp;├── orbitsd{d}_n{n}_separated_sorted/ <- These folders include all orbits for d the local dimension and n the number of particles. <br>
     &emsp;&emsp;&emsp;&emsp;├── orbit_i <- The index i indicates the different orbits.         
&emsp;&emsp;├── aspl.ipynb <- Notebook for the extraction of the average shortest path and the heuristic method adopted. <br>
&emsp;&emsp;├── auxilary_functions.py <- Helper functions for qudit_lc.py <br>
&emsp;&emsp;├── data_analysis_functions.py <- A module with every function required for the data extraction and analysis. <br>
&emsp;&emsp;├── main_data_notebook.ipynb <- The notebook with the statistical analysis after the data extraction. <br>
&emsp;&emsp;├── min_edge_color_states_wolfram_engine.py <- Module to compute the minimum edge chromatic number. <br>
&emsp;&emsp;├── mspl.ipynb <- Notebook for the extraction of the maximum shortest path and the heuristic method adopted. <br>
&emsp;&emsp;├── qudit_lc.py <- Main functions for graph generation, orbit extraction, and Schmidt measure computation. <br>
&emsp;&emsp;├── real_join_data.py <- Module for the extraction of the required datasets. <br>
&emsp;&emsp;└── sm_final_data.py <- Module for preprocessing the Schmidt Measure data before the data analysis. <br>

