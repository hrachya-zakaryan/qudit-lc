"""
Script to compute the minimum edge chromatic number (edge coloring) 
of graph orbits represented as adjacency matrices.

This script loads graph data files from Mathematica-generated directories,
parses them using a Wolfram Language session, constructs graphs, and computes
the edge chromatic number for each. The minimum edge chromatic number for
each orbit is collected and written to a summary file.

Dependencies:
- networkx
- numpy
- os
- WolframClient for Python (wolframclient)
- Custom modules: data_analysis_functions, isomorphism_method

Inputs:
- Directory structure: "files_for_mathematica/orbits_d{d}_n{cnt}_separated"
  Each file should contain a list of adjacency matrices in Mathematica syntax.

Outputs:
- A text file "og_data/join_min_edge_colors_multigraphs_in_OG.txt" containing
  a Python list of the minimum edge chromatic numbers for all orbits.

Usage:
- Modify `d` and `n` to change the graph dimension and size.
- Ensure the Wolfram Engine is installed and accessible.

Notes:
-The output is not used in the paper because the maximum degree is used. The reason for this was mainly that we wanted all the observables to be calculated using Python and software that is not accessible unless one has a license. On the other hand, if the user has a Wolfram Engine can use our short script directly.
"""


import networkx as nx
import numpy as np
from data_analysis_functions import *
from isomorphism_method import *
import os
import networkx as nx
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl

# Start a Wolfram Engine session
session = WolframLanguageSession()

d = 3
n = 7

folder_base = "files_for_mathematica"
final_array_of_min_edge_color = []

for cnt in range(3, n + 1):
    print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
    
    # Specify the correct directory path
    orbit_folder = f"{folder_base}/orbits_d{d}_n{cnt}_separated"
    
    # Get the list of files in the orbit folder
    files_for_given_orbit = os.listdir(orbit_folder)
    
    # Sort the files based on the orbit number extracted from the filenames
    files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1].replace('.txt', '')))
    
    for file_path in files_for_given_orbit:
        # Generate the full file path by combining the folder and file name
        full_file_path = os.path.join(orbit_folder, file_path)
        
        # Import the text file as lines using the full path
        raw_text = session.evaluate(wl.Import(full_file_path, "Lines"))
        
        # Remove empty entries (Use Wolfram's built-in function)
        filtered_text = session.evaluate(wl.DeleteCases(raw_text, ""))
        
        # Convert text to Mathematica expressions (matrices)
        matrices = session.evaluate(wl.Map(wl.ToExpression, filtered_text))
        
        color_of_edges_for_specific_orbit = []
        for r in matrices:
            # Graph processing
            calculate_graph = session.evaluate(wl.AdjacencyGraph(r))
            edge_chromatic_number = session.evaluate(wl.EdgeChromaticNumber(calculate_graph))
            color_of_edges_for_specific_orbit.append(edge_chromatic_number)
        
        # Get the minimum edge color for this orbit
        min_edge_color_orbit_specific = min(color_of_edges_for_specific_orbit)
        final_array_of_min_edge_color.append(min_edge_color_orbit_specific)
        print(f"Min edge color for {file_path}: {min_edge_color_orbit_specific}")
        
# Close the session
session.terminate()

#write an external file
with open("og_data/join_min_edge_colors_multigraphs_in_OG.txt", "w") as f:
    f.write(str(final_array_of_min_edge_color))
