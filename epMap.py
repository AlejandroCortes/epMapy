###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
import numpy as np  #do maths on multidimensional array objects
import pandas as pd  #to read data sets
import re  #to compare a set of string-type attributes
import matplotlib.pyplot as plt #default plotting library
# importing external functions from functionsmapping.py
from functionsepmapy import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data, extract_profile, show_intro  
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
import scipy.ndimage
import os # to access paths and modify file names
import json # to work with json files
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
data = load_data(file_path, file_type, data_structure) # from the variables extracted from the chosen file, loads the data
print("data file was read") # to keep track on things

###########################################################################################################################################################
#Structuring data so that it can be plotted
###########################################################################################################################################################
def create_oxide_grids(data): #function to build grids of oxides for easy plotting based on NX and NY markers
    oxide_grids = {} #empty dictionary

    for key, value in data.items(): # Initialise the grid for each oxide present in the data
        oxides_data = extract_elementsoxides(value) # Extract oxides for each row using the existing function in FileSelector
        for oxide, oxide_value in oxides_data.items(): # For each oxide present in the row, create an entry in the grid if not already present
            if oxide not in oxide_grids:
                oxide_grids[oxide] = {}
            NX = value['NX'] # get the NX marker of the pixel
            NY = value['NY'] # get the NY marker of the pixel
            X_coord = value.get('X_coord', None)  # get real X coordinates in mm
            Y_coord = value.get('Y_coord', None)  # get real Y coordinates in mm
            if (NX, NY) not in oxide_grids[oxide]: # Initialise the grid based on marker combination (NX, NY) for the oxide
                oxide_grids[oxide][(NX, NY)] = {'oxide_value': None, 'X_coord': None, 'Y_coord': None}  # Set to None initially
            oxide_grids[oxide][(NX, NY)] = {'oxide_value': oxide_value, 'X_coord': X_coord, 'Y_coord': Y_coord} #filling values

    return oxide_grids

oxide_grids = create_oxide_grids(data) # returns a dictionary in grid like distribution
new_file_name = re.sub(file_type, '', os.path.basename(file_path))
sample_name = new_file_name.strip()
def convert_oxide_grids_for_json(oxide_grids): #function to save the oxide_grids as a json file
    new_oxide_grids = {} # empty dictionary to overcome the issue of having two markers (NX,NY) instead of one
    for oxide, grid in oxide_grids.items():
        new_grid = {} #empty subgrid per oxide
        for coordinate, value in grid.items(): #puts together Nx and Ny with and underscore to solve issues
            coordinate_str = f"{coordinate[0]}_{coordinate[1]}"  # Convert (NX, NY) to 'NX_NY'
            #store the oxide value along with its coordinates
            new_grid[coordinate_str] = {
                'oxide_value': value['oxide_value'],
                'X_coord': value['X_coord'],
                'Y_coord': value['Y_coord']
            }

        new_oxide_grids[oxide] = new_grid #add the grid to the main oxide_grids
    with open(sample_name+'oxide_grids.json', 'w') as f:
        json.dump(new_oxide_grids, f, indent=4)
    print("Oxide grids saved to JSON file.")


oxide_grids_str_keys = convert_oxide_grids_for_json(oxide_grids) # saves the json file




# # Initialise arrays for various properties
# # molar indexes: NK_A:Na2O+K2O/Al2O3, NKC_A:Na2O+K2O+CaO/Al2O3, K_Na:K2O/Na2O, M_f:
# NK_A, NKC_A, Mg_MgFe2, K_Na, M_f = [np.zeros((c, d)) for _ in range(5)]

# #creating grids with coordinates and
# coord_grids = {coordinate: np.zeros((c,d)) for coordinate in coord_column_indices.keys()}
# oxide_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_original.keys()}
# oxidef_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_original.keys()}
# data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
# print(oxide_column_indices)
# fe2o3_included, oxide_column_indices_fe = add_fe2o3(data_majors,oxide_column_indices) #calculates fe2o3 and updates the oxide_column_indices
# print(oxide_column_indices_fe,oxide_column_indices)
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
