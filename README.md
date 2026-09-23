# qudit-lc

A complete framework for calculating the Local Clifford (LC) entanglement classes for qudit graph states, the necessary observables extraction, and the corresponding data analysis.

## Overview

The goal of this work was the extraction of the entanglement classes using local complementation.

- **Orbit extraction:** For this step, we had to calculate the LC orbits and represent them as a graph called Orbit Graph (OG).
- **Define Observables:** At this step, we defined a set of properties emerging from graph theory as well as from quantum theory. For example, the vertex number, the density of the graph, and the Schmidt Measure.
- **Correlation Analysis and Results Presentation:** After the extraction of the datasets, we calculated the Pearson and Kendall's tau coefficients to determine if some of the obtained datasets are correlated. Then some 
Scripts are provided for presenting our results with figures and tables for latex.
- **Benchmarking:** Our code is cross-checked with results already known in the literature for qubits. The results are in agreement with Quantum 4, 305 (2020).

## Installation

```git clone https://github.com/hrachya-zakaryan/qudit-lc```


## Python Version and Packages
The user must ensure to use the following versions of the packages.<br>
Python                    3.11.2 <br>
matplotlib                3.11.1 <br>
numpy                     2.4.6 <br>
pandas                    3.0.5 <br>
scipy                     1.17.1 <br>
ipykernel                 7.3.0 (only needed to run the notebook) <br>

The two C programs need a C compiler with OpenMP support (tested with gcc 12.2.0). orbits_mt is built with the line in its header,<br>
```cc -O3 -march=native -fopenmp -o orbits_mt orbits_mt.c```<br>
and orbit_scan is compiled automatically by scan_tool.py on first use.

## File Hierarchy
qudit-lc/ <br>
&emsp;&emsp;├── orbits_mt.c <- Parallel LC-orbit enumeration for qudit graph states (OpenMP, lock-free union-find). <br>
&emsp;&emsp;├── orbit_scan.c <- One streaming pass per level: orbit representative, member observables and Schmidt-measure bounds. <br>
&emsp;&emsp;├── scan_tool.py <- Builds and runs orbit_scan and parses its output. <br>
&emsp;&emsp;├── fast_analysis.py <- Orbit-graph observables on a scipy CSR adjacency: ASPL, diameter, density, chromatic bound. <br>
&emsp;&emsp;├── make_table.py <- Per-orbit table for every level: parallel driver, csv and LaTeX longtable output. <br>
&emsp;&emsp;├── representatives.py <- Extraction and drawing of the orbit representatives. <br>
&emsp;&emsp;├── og_analysis_fast.ipynb <- Notebook driving the pipeline and the data analysis. <br>
&emsp;&emsp;├── data_analysis.ipynb <- Notebook performing the statistical analysis. <br>
&emsp;&emsp;├── table_d2_n9.csv <- Per-orbit data for d = 2, n = 3..9 (585 orbits). <br>
&emsp;&emsp;├── table_d3_n8.csv <- Per-orbit data for d = 3, n = 3..8 (762 orbits). <br>
&emsp;&emsp;└── table_d5_n6.csv <- Per-orbit data for d = 5, n = 3..6 (49 orbits). <br>
