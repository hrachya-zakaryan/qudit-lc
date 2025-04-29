"""
Script to flatten and modify nested lists by either calculating the mean of sublists with 
two elements or converting them into tuples.

This script defines functions to recursively flatten nested lists and modify sublists with
exactly two elements. The `flatten_and_modify` function computes the mean of these sublists, 
while the `flatten_and_modify_table` function replaces them with tuples. The script also
includes a `main` function to load data from a file, flatten it using one of the methods, 
and print the resulting list.

Both of the functions are useful for our work. The first one calculates the array used in the
data analysis. The second one creates the appropriate array to produce the table 
for the values of each observable for every orbit. This is relevant to the Schmidt Measure.

Dependencies:
- numpy

Inputs:
- A text file "og_data/join_sm_in_og.txt" containing a nested list in Python literal format.

Outputs:
- A flattened list printed to the console, with sublists of length 2 replaced by either
  their mean (for `flatten_and_modify`) or a tuple (for `flatten_and_modify_table`).

Usage:
- Modify the `RUN_MAIN` variable to `True` to run the loading, flattening, and printing steps.

Notes:
- Ensure that the file path "og_data/join_sm_in_og.txt" is correct and contains a valid nested list.
- The `RUN_MAIN` flag allows the user to control when the script should execute the main actions.
"""
import numpy as np

def flatten_and_modify(nested_list):
    """
    Flattens a nested list, replacing sublists with exactly two elements by their mean.

    This function recursively flattens a nested list of arbitrary depth. 
    When encountering a sublist of length 2, it replaces the sublist with the mean 
    of its two elements. All other elements are added to the flattened list unchanged.

    Parameters:
    nested_list (list): The nested list to be flattened and modified.

    Returns:
    list: A new list that is flattened and has sublists of length 2 replaced by their mean.

    Example:
    >>> flatten_and_modify([[1, 2], [3, 4], 5])
    [1.5, 3.5, 5]
    """
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):  # If item is a list, check if it needs modification
            if len(item) == 2:  # If the list has exactly two elements, take their mean
                flat_list.append(np.mean(item))
            else:
                flat_list.extend(flatten_and_modify(item))  # Recursively flatten and add to the list
        else:
            flat_list.append(item)  # If it's not a list, just append it
    return flat_list


def flatten_and_modify_table(nested_list):
    """
    Flattens a nested list, replacing sublists with exactly two elements by a tuple (a, b).

    This function recursively flattens a nested list of arbitrary depth. 
    When encountering a sublist of length 2, it replaces the sublist with a tuple 
    containing the two elements. All other elements are added to the flattened list unchanged.

    Parameters:
    nested_list (list): The nested list to be flattened and modified.

    Returns:
    list: A new list that is flattened and has sublists of length 2 replaced by a tuple.

    Example:
    >>> flatten_and_modify_table([[1, 2], [3, 4], 5])
    [(1, 2), (3, 4), 5]
    """
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            if len(item) == 2:
                flat_list.append(tuple(item))  # Replace sublist with tuple
            else:
                flat_list.extend(flatten_and_modify_table(item))  # Recursively flatten
        else:
            flat_list.append(item)  # If it's not a list, just append it
    return flat_list


# Set this variable to True if one wants to run the code below
RUN_MAIN = False

def main():
    """
    Main function to load data, flatten it using the appropriate method, and print the result.

    This function opens a file containing a nested list, flattens the list using the
    `flatten_and_modify_table` function, and prints the resulting flattened list.

    Example:
    >>> main()  # This will load the data, flatten it, and print the result
    """
    # Load the data from the file
    with open("og_data/join_sm_in_og.txt", "r") as f:
        data_sm_initial = eval(f.read())
    
    # Flatten the data with modifications
    flattened_data = flatten_and_modify_table(data_sm_initial)
    
    print(flattened_data)


if __name__ == "__main__":
    if RUN_MAIN:
        main()
