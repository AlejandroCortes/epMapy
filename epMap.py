###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
import numpy as np  #do maths on multidimensional array objects
import pandas as pd  #to read data sets
import re  #to compare a set of string-type attributes
import matplotlib.pyplot as plt #default plotting library
# importing external functions from functionsmapping.py
from functionsepmapy import to_anhydrous, k_na, to_mol, to_cat, norm_calc, add_fe2o3, clean_data, show_intro, load_json, m_f  
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
import scipy.ndimage
import os # to access paths and modify file names
import json # to work with json files
import pickle
# functions and stuff from FileSelector.py
from file_selector import FileSelector, load_data, extract_elementsoxides 

###########################################################################################################################################################
#  Intro messages
###########################################################################################################################################################
welcome_message = """
Input file with major and minor element composition"""
instructions_message = """
text (tab-separated) or xlsx file with mapping grid
"""
InstructionWindow = FileSelector(welcome_message, instructions_message) #creates window and gives the option to browse a file
file_path, file_type = InstructionWindow.get_file_path_and_type() #obtains the file type and path from the browse window

###########################################################################################################################################################
# Loading EPMA data
###########################################################################################################################################################
data_structure = 'map' # this is for the FileSelector to know how to read the file (map or single-point)
print("processing "+data_structure+" data stored in "+os.path.basename(file_path)+" file...")
data = load_data(file_path, file_type, data_structure) # from the variables extracted from the chosen file, loads the data
print("data file was read") # to keep track on things

###########################################################################################################################################################
#Structuring data so that it can be plotted
###########################################################################################################################################################
def create_oxide_grids(data):
    oxide_grids = {}  # Empty dictionary to store oxide grids
    total_grids = {}  # Separate dictionary to store the 'Total' oxide values

    for _, value in data.items():
        oxides_data = extract_elementsoxides(value)  # Extract oxides for each sample
        total = value.get('Total', 0)  # Assuming 'TOTAL' is the column in the data that has the total oxide weight

        for oxide, oxide_value in oxides_data.items():  # Iterate through extracted oxides
            if oxide not in oxide_grids:
                oxide_grids[oxide] = {}

            NX = value['NX']  # Get NX value
            NY = value['NY']  # Get NY value
            X_coord = value.get('X_coord', None)  # Get X coordinate (if available)
            Y_coord = value.get('Y_coord', None)  # Get Y coordinate (if available)

            # If this is the first time we are encountering this (NX, NY) pair for the current oxide
            if (NX, NY) not in oxide_grids[oxide]:
                oxide_grids[oxide][(NX, NY)] = {
                    'oxide_value': None,
                    'X_coord': X_coord,
                    'Y_coord': Y_coord
                }

            oxide_grids[oxide][(NX, NY)] = {
                'oxide_value': oxide_value,
                'X_coord': X_coord,
                'Y_coord': Y_coord
            }

        # Add Total as a separate "oxide" entry in the total_grids dictionary
        if (NX, NY) not in total_grids:
            total_grids[(NX, NY)] = {}

        total_grids[(NX, NY)] = total  # Store the Total value for this (NX, NY) pair

    return oxide_grids, total_grids

oxide_grids, total_grids = create_oxide_grids(data) # returns a dictionary in grid like distribution
clean_grids, total_grids = clean_data(oxide_grids,total_grids)
anhydrous_grids, total_grids = to_anhydrous(clean_grids, total_grids)
k_na_grid = k_na(anhydrous_grids, total_grids)
new_file_name = re.sub(file_type, '', os.path.basename(file_path))
sample_name = new_file_name.strip()

# # options = [
# #     #("Plot K/Na", function_1),
# #     #("Plot NaK/Al", function_2),
# #     #("Plot NaKCa/Al", function_3),
# #     #("Plot zrc saturation", function_4),
# #     #("Plot normative mineralogy", function_5),
# #     ("Plot raw oxides", clean_data(data)),
# #     #("Plot anhydrous oxides", function_7),
# #     ("save raw oxides as grids in json format", convert_oxide_grids_for_json)
# # ]
# # oxide_grids_str_keys = convert_oxide_grids_for_json(oxide_grids) # saves the json file

# #Initialise arrays for various properties
# # molar indexes: NK_A:Na2O+K2O/Al2O3, NKC_A:Na2O+K2O+CaO/Al2O3, K_Na:K2O/Na2O, M_f:

# data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
# fe2o3_included, oxide_column_indices_fe = add_fe2o3(data_majors,oxide_column_indices) #calculates fe2o3 and updates the oxide_column_indices
# oxidef_grids_fe = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_fe.keys()}
# data_anhf = to_anhydrous(data_majors, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
# data_anhf_fe = to_anhydrous(fe2o3_included, oxide_column_indices_fe) #calculates anhydrous-base oxide compositions wt.%
# data_molf   = to_mol(data_anhf, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
# data_molf_fe   = to_mol(data_anhf_fe, oxide_column_indices_fe)         #calculates anhydrous-based oxide compositions mol
# data_catf   = to_cat(data_molf, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
# data_mol   = to_mol(data_majors, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
# data_cat   = to_cat(data_mol, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
# data_normf  = norm_calc(data_molf_fe, oxide_column_indices_fe) #calculates normative mineralogy

# #oxidefe3_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}

# #molanh_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}

# # Loop through all data points to fill grid arrays
# for idx in range(len(coord_data[coord_column_indices["NXY"]])):
#     nx, ny = int(coord_data[coord_column_indices["NX"]][idx]) - 1, int(coord_data[coord_column_indices["NY"]][idx]) - 1
    
#     if data_majors[oxide_column_indices["Total"]][idx] > 0:
#         for oxide, col_index in oxide_column_indices_original.items():
#             col_index = oxide_column_indices[oxide]
#             oxide_grids[oxide][nx, ny] = data_majors[col_index][idx]
#             oxidef_grids[oxide][nx, ny] = data_anhf[col_index][idx]
#             oxidef_grids_fe[oxide][nx, ny] = data_anhf_fe[col_index][idx]
#         for coordinate, col_index in coord_column_indices.items():
#             col_index = coord_column_indices[coordinate]
#             coord_grids[coordinate][nx, ny] = coord_data[col_index][idx]

#         # Fill derived properties
#         M_f[nx, ny] = (
#             (data_catf[oxide_column_indices["Na2O"],idx] +
#              data_catf[oxide_column_indices["K2O"],idx] +
#              (data_catf[oxide_column_indices["CaO"],idx] * 2)) /
#             (data_catf[oxide_column_indices["SiO2"],idx] *
#              data_catf[oxide_column_indices["Al2O3"],idx])
#         )
#         NK_A[nx, ny] = (
#             (data_molf[oxide_column_indices["Na2O"],idx] +
#              data_molf[oxide_column_indices["K2O"],idx]) /
#             data_molf[oxide_column_indices["Al2O3"],idx]
#         )
#         NKC_A[nx, ny] = (
#             (data_molf[oxide_column_indices["Na2O"],idx] +
#              data_molf[oxide_column_indices["K2O"],idx] +
#              data_molf[oxide_column_indices["CaO"],idx]) /
#             data_molf[oxide_column_indices["Al2O3"],idx]
#         )
#         K_Na[nx, ny] = (
#             data_anhf[oxide_column_indices["K2O"],idx] * 8302 /
#             data_anhf[oxide_column_indices["Na2O"],idx] * 7419
#         )
#     else:
#         for oxide in oxide_column_indices_original.keys():
#             oxide_grids[oxide][nx, ny] = -1  # Placeholder for invalid points
#             oxidef_grids[oxide][nx, ny] = -1  # Placeholder for invalid points
#         for grid in [M_f, NK_A, NKC_A, K_Na]:
#             grid[nx, ny] = -1
#         for coordinate, col_index in coord_column_indices.items():
#             col_index = coord_column_indices[coordinate]
#             coord_grids[coordinate][nx, ny] = coord_data[col_index][idx]


# ###########################################################################################################################################################
# # Plotting and profile extraction loop
# ###########################################################################################################################################################
# # Main loop for plotting and profile extraction

# pixel_size = (coord_data[coord_column_indices["X"], 1] - coord_data[coord_column_indices["X"], 0]) * 1000

# while True:
#     selected_oxides_input = input("Enter the oxides you want to plot, separated by commas, or type 'all' to plot all: ").strip()
#     if selected_oxides_input.lower() == 'all':
#         selected_oxides = list(oxide_columns_original.keys())  # Select all available oxides
#     else:
#         selected_oxides = [oxide.strip() for oxide in selected_oxides_input.split(',')]  # Creates a list with the selected oxides

#     cm = 1 / 2.54
#     plt.figure(figsize=(6.5 * cm, 6.5 * cm))
#     fig, axs = plt.subplots(1, len(selected_oxides), figsize=(3 * len(selected_oxides), 3))

#     for i, oxide in enumerate(selected_oxides):
#         oxide = oxide.strip()
#         if oxide in oxide_grids:
#             ax = axs[i] if len(selected_oxides) > 1 else axs
#             im = ax.imshow(oxide_grids[oxide], cmap='viridis', vmin=0, vmax=np.max(oxide_grids[oxide]))
#             im.cmap.set_under('black')
#             ax.set_title(f"{oxide} wt.%")
#             ax.axis("off")
#             fig.colorbar(im, ax=ax)
#         else:
#             print(f"Warning: {oxide} is not a valid oxide.")

#     scalebar = ScaleBar(pixel_size, "um")
#     plt.gca().add_artist(scalebar)
#     plt.tight_layout()

#     formatted_oxides = "_".join([oxide.strip() for oxide in selected_oxides])
#     pdf_filename = f"{sample_name}_{formatted_oxides}.pdf"
#     plt.savefig(pdf_filename, dpi=600, transparent=True, bbox_inches='tight')
#     print(f"Plot saved as {pdf_filename}")
#     plt.show()

#     # Now, enter a loop to allow multiple profile extractions
#     profile_index = 1  # Initialize profile index
#     while True:
#         # Ask User if they want to extract a profile
#         extract = input("Would you like to extract a traverse for the selected oxides? (yes/no): ").strip().lower()
#         if extract == 'yes':
#             # Ask user which oxide to use for selecting the profile points
#             selected_oxide = input(f"Which oxide would you like to use for selecting the profile points? Choose from: {', '.join(selected_oxides)}: ").strip()

#             if selected_oxide not in selected_oxides:
#                 print("Invalid oxide selection. Please select an oxide from the list.")
#                 continue

#             # Extract profile for the selected oxide
#             extract_profile(oxide_grids, selected_oxides, pixel_size, selected_oxide, sample_name, profile_index)
#             profile_index += 1  # Increment the profile index after each extraction

#         else:
#             break  # Exit the profile extraction loop if the user is done

#     # Ask if the user wants to plot another set of oxides after extracting profiles
#     repeat = input("Would you like to plot another set of oxides? (yes/no): ").strip().lower()
#     if repeat != 'yes':
#         break  # Exit the entire loop if the user doesn't want to plot another set

def plot_oxide_grids(oxide_grids, sample_name):
    selected_oxides = list(oxide_grids.keys())  # Automatically select all oxides in the oxide_grids dictionary

    ## Get the first two consecutive rows to calculate the pixel size (X_coord, Y_coord)
    first_oxide = selected_oxides[0]  # Take the first oxide (assuming all oxides have the same grid)
    first_data = list(oxide_grids[first_oxide].items())  # Get the data for the first oxide

    # Take the first two consecutive pixels (using X_coord, Y_coord) for the first oxide
    (NX1, NY1), data1 = first_data[0]
    (NX2, NY2), data2 = first_data[1]

    # Get the X_coord and Y_coord for the first two data points
    x1, y1 = data1['X_coord'], data1['Y_coord']
    x2, y2 = data2['X_coord'], data2['Y_coord']

    # Calculate pixel size by the distance between the first two pixels in the X or Y direction
    pixel_size_x = abs(x2 - x1)  # Calculate pixel size based on X_coord
    pixel_size_y = abs(y2 - y1)  # Calculate pixel size based on Y_coord
    pixel_size = np.sqrt(pixel_size_x**2 + pixel_size_y**2)*1000  # Overall pixel size, assuming the grid is square
    
    # Calculate the number of rows and columns for the subplot layout
    num_oxides = len(selected_oxides)
    num_columns = 3
    num_rows = (num_oxides + num_columns - 1) // num_columns  # Round up to ensure enough rows

    fig, axs = plt.subplots(num_rows, num_columns, figsize=(3 * num_columns, 3 * num_rows))
    axs = axs.flatten()  # Flatten the axes array to make indexing easier

    for i, oxide in enumerate(selected_oxides):
        oxide = oxide.strip()
        if oxide in oxide_grids:
            ax = axs[i]  # Select the correct axis for the oxide plot

            max_NX = max([key[0] for key in oxide_grids[oxide].keys()])  # Create a grid based on NX, NY values
            max_NY = max([key[1] for key in oxide_grids[oxide].keys()])
            oxide_matrix = np.full((max_NX + 1, max_NY + 1), np.nan)  # Initialize an empty matrix with NaN values for plotting
            
            for (NX, NY), data in oxide_grids[oxide].items():  # Populate the matrix with oxide values
                oxide_matrix[NX, NY] = data['oxide_value']
                
            im = ax.imshow(oxide_matrix, cmap='viridis', vmin=0, vmax=np.nanmax(oxide_matrix))  # Plot the oxide matrix with imshow and apply colormap
            im.cmap.set_under('black')  # Set color for NaN values
            ax.set_title(f"{oxide} wt.%")
            ax.axis("off")  # Hide axis
            fig.colorbar(im, ax=ax)
            scalebar = ScaleBar(pixel_size, "um", location="lower right")
            ax.add_artist(scalebar)
 
        else:
            print(f"Warning: {oxide} is not a valid oxide.")
    
    # Hide any unused subplots
    for j in range(i + 1, len(axs)):
        axs[j].axis('off')

    # Add scale bar
    
    
    
    plt.tight_layout()

    formatted_oxides = "_".join([oxide.strip() for oxide in selected_oxides])
    extraname = r'_FirstPass_\d{5}__Oxide_Image_Classify'  # Part of the file name to remove
    new_sample_name = re.sub(extraname, '', sample_name)  # Simplify the file name
    pdf_filename = f"{new_sample_name}_{formatted_oxides}.pdf"

    plt.savefig(pdf_filename, dpi=600, transparent=True, bbox_inches='tight')
    print(f"Plot saved as {pdf_filename}")
    plt.show()
    return pixel_size, new_sample_name



pixel_size, new_sample_name = plot_oxide_grids(clean_grids, sample_name)


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

# Global variables to store selected points
selected_points = []

def on_click(event, oxide_grids, oxide, pixel_size):
    """Callback function for mouse click events"""
    global selected_points

    # Check if the click is within the image
    if event.inaxes is not None:
        # Store the coordinates of the click (NX, NY)
        x = int(event.xdata)
        y = int(event.ydata)

        # Append the selected point
        selected_points.append((x, y))

        # Plot the selected point on the image
        event.inaxes.plot(x, y, 'ro')  # Mark the point with a red 'o'
 
        if len(selected_points) == 2:
            # Once two points are selected, extract and plot the profile
            extract_and_plot_profile(oxide_grids, oxide, selected_points[0], selected_points[1], pixel_size)

        plt.draw()  # Redraw the plot to update it

def extract_and_plot_profile(oxide_grids, oxide, point1, point2, pixel_size):
    """Extract and plot the oxide profile between two points in real distance (micrometers)"""
    # Get the coordinates for the two points
    x1, y1 = point1
    x2, y2 = point2

    # Get the real-world coordinates for the selected points (convert mm to um by multiplying by 1000)
    coord1 = oxide_grids[oxide].get((x1, y1), {})
    coord2 = oxide_grids[oxide].get((x2, y2), {})

    if not coord1 or not coord2:
        print("Error: One or both points are not in the oxide grid.")
        return

    # Extract the X and Y coordinates for both points (convert from mm to µm)
    x1_real, y1_real = coord1['X_coord'] * 1000, coord1['Y_coord'] * 1000  # Convert to micrometers
    x2_real, y2_real = coord2['X_coord'] * 1000, coord2['Y_coord'] * 1000  # Convert to micrometers

    print(f"Selected Points: ({x1_real}, {y1_real}) and ({x2_real}, {y2_real})")  # Debug print

    # Calculate the Euclidean distance between the points in real-world space (micrometers)
    dist = np.sqrt((x2_real - x1_real)**2 + (y2_real - y1_real)**2)
    print(f"Real-world distance: {dist} µm")  # Debug print

    # Calculate the number of points based on the pixel size
    num_points = int(dist / pixel_size)
    
    # Limit the number of points if it's too high (for performance reasons)
    if num_points > 500:
        print("Too many points along the profile, limiting to 500 points.")
        num_points = 500

    print(f"Number of points along the profile: {num_points}")  # Debug print

    # Generate points along the line between the two selected points
    profile_x = np.linspace(x1_real, x2_real, num_points)
    profile_y = np.linspace(y1_real, y2_real, num_points)

    # Prepare a list of the (NX, NY) coordinates and their corresponding oxide values
    grid_points = []
    oxide_values = []
    
    for (NX, NY), data in oxide_grids[oxide].items():
        # Store the grid points and oxide values
        grid_points.append((data['X_coord'] * 1000, data['Y_coord'] * 1000, data['oxide_value']))  # Convert to um
    grid_points = np.array(grid_points)

    # Use KDTree for efficient nearest-neighbor search
    kdtree = KDTree(grid_points[:, :2])  # Only use the X and Y coordinates for search

    profile_values = []
    real_distances = []  # To store real-world distances

    # Iterate over the profile line points and get the oxide values
    for px, py in zip(profile_x, profile_y):
        # Find the closest (NX, NY) in the grid using KDTree for faster search
        dist, index = kdtree.query([px, py])  # Get the nearest point in the grid
        oxide_value = grid_points[index, 2]  # Get the oxide value at the closest grid point
        profile_values.append(oxide_value)

        # Compute the real-world distance from the first point (in µm)
        real_distance = np.sqrt((px - x1_real)**2 + (py - y1_real)**2)
        real_distances.append(real_distance)

    print(f"Profile values: {profile_values}")  # Debug print
    print(f"Real distances: {real_distances}")  # Debug print

    # Plot the oxide profile against real-world distance (in micrometers)
    plt.figure()
    plt.plot(real_distances, profile_values, label=f'{oxide} profile')
    plt.title(f"Profile of {oxide} between points ({x1_real}, {y1_real}) and ({x2_real}, {y2_real})")
    plt.xlabel('Distance along the profile (µm)')
    plt.ylabel(f'{oxide} wt.%')
    plt.legend()
    plt.show()

def plot_oxide_grids_interactive(oxide_grids, sample_name, pixel_size):
    """Plot the oxide grids and allow the user to select two points for a profile"""
    selected_oxides = list(oxide_grids.keys())

    # Ask the user to select an oxide
    print("Available oxides: ", selected_oxides)
    oxide = input("Please enter the oxide you would like to plot the profile for: ").strip()

    # Check if the entered oxide is valid
    if oxide not in selected_oxides:
        print(f"Error: '{oxide}' is not a valid oxide. Please choose one from the available oxides.")
        return  # Exit if the oxide is not valid

    # Get the grid dimensions for the selected oxide
    max_NX = max([key[0] for key in oxide_grids[oxide].keys()])
    max_NY = max([key[1] for key in oxide_grids[oxide].keys()])
    oxide_matrix = np.full((max_NX + 1, max_NY + 1), np.nan)

    # Populate the matrix with oxide values
    for (NX, NY), data in oxide_grids[oxide].items():
        oxide_matrix[NX, NY] = data['oxide_value']

    # Plot the oxide matrix
    fig, ax = plt.subplots()
    im = ax.imshow(oxide_matrix, cmap='viridis', vmin=0, vmax=np.nanmax(oxide_matrix))
    ax.set_title(f"{oxide} wt.%")
    ax.axis("off")
    fig.colorbar(im)

    # Add a cursor to show the mouse coordinates
    cursor = Cursor(ax, useblit=True, color='red', linewidth=1)

    # Connect the click event to select points
    cid = fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, oxide_grids, oxide, pixel_size))

    plt.show()

# Usage:
# Assuming oxide_grids is already created
# sample_name = "your_sample_name"
# pixel_size = 0.5  # Example pixel size in micrometers
#plot_oxide_grids_interactive(clean_grids, new_sample_name, pixel_size)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import re

def plot_k_na_grid(k_na_grid, mf_grid, sample_name):
    # Extract unique NX, NY coordinates
    coordinates = list(k_na_grid.keys())
    max_NX = max([coord[0] for coord in coordinates])
    max_NY = max([coord[1] for coord in coordinates])

    # Create a matrix to hold the k_na values (with NaN for missing data)
    k_na_matrix = np.full((max_NX + 1, max_NY + 1), np.nan)
    mf_matrix = np.full((max_NX + 1, max_NY + 1), np.nan)

    # Populate the k_na matrix with the k_na values
    for (NX, NY), data in k_na_grid.items():
        k_na_matrix[NX, NY] = data['oxide_value']
    for (NX, NY), data2 in mf_grid.items():
        mf_matrix[NX, NY] = data2['oxide_value']
    

    # Calculate pixel size based on the X_coord and Y_coord of the first two data points
    first_data = list(k_na_grid.values())[0]
    x1, y1 = first_data['X_coord'], first_data['Y_coord']
    x2, y2 = list(k_na_grid.values())[1]['X_coord'], list(k_na_grid.values())[1]['Y_coord']

    pixel_size_x = abs(x2 - x1)  # Calculate pixel size based on X_coord
    pixel_size_y = abs(y2 - y1)  # Calculate pixel size based on Y_coord
    pixel_size = np.sqrt(pixel_size_x**2 + pixel_size_y**2) * 1000  # Overall pixel size in microns (assuming grid is square)

    # Create the plot
    fig, axs = plt.subplots(1,2,figsize=(6, 3))

    ax = axs[0]

    # Use imshow to display the k_na matrix with the 'viridis' colormap
    im = ax.imshow(k_na_matrix, cmap='viridis', vmin=0, vmax=6)  # Set range from 0 to max k_na value
    im.cmap.set_under('black')  # Set color for NaN values (black)
    
    # Add title and colorbar
    ax.set_title("K/Na molar ratio")
    ax.axis('off')  # Hide axes
    fig.colorbar(im, ax=ax)
    
    # Add scale bar
    scalebar = ScaleBar(pixel_size, "um", location="lower right")
    ax.add_artist(scalebar)

    ax = axs[1]

    # Use imshow to display the k_na matrix with the 'viridis' colormap
    im = ax.imshow(mf_matrix, cmap='viridis', vmin=1, vmax=2)  # Set range from 0 to max k_na value
    im.cmap.set_under('black')  # Set color for NaN values (black)
    
    # Add title and colorbar
    ax.set_title("M cation ratio")
    ax.axis('off')  # Hide axes
    fig.colorbar(im, ax=ax)
    
    # Add scale bar
    scalebar = ScaleBar(pixel_size, "um", location="lower right")
    ax.add_artist(scalebar)

    # Save the plot as a PDF file
    formatted_name = sample_name.replace('FirstPass', '').replace(r'\d{5}', '')  # Clean up the sample name
    pdf_filename = f"{formatted_name}_KNa_Mf_grid_plot.pdf"
    plt.savefig(pdf_filename, dpi=600, transparent=True, bbox_inches='tight')

    # Show the plot
    plt.show()

    print(f"Plot saved as {pdf_filename}")
    return pixel_size, formatted_name




elements_weights, oxides_weights = load_json()

mol_grids, total_grids = to_mol(anhydrous_grids, total_grids, oxides_weights)

cat_grids, total_grids = to_cat(mol_grids, total_grids, oxides_weights)

mf_grid, total_grids = m_f(cat_grids, total_grids)

plot_k_na_grid(k_na_grid, mf_grid, new_sample_name)