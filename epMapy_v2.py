###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
import numpy as np  #do maths on multidimensional array objects
import pandas as pd  #to read data sets
import re  #to compare a set of string-type attributes
import time  #to handle time-related calculations
import matplotlib.pyplot as plt #default plotting library
# importing external functions from functionsmapping.py
from functionsepmapy import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data, extract_profile, show_intro  
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
import scipy.ndimage
import os

###########################################################################################################################################################
# Important information for User
###########################################################################################################################################################
Introd = """
epMapy allows plotting EPMA quantitative maps and performing calculations with the data within, 
as well as treatment of single-point data acquired with the microprobe.

Please be aware that this script is designed for the output file obtained with CalcImage (J. Donovan) when processing maps.

The input file for maps and single-point analyses should be in .DAT format 

Regarding maps, the heading of the file should contain at least the following column names (character sensitive): X, Y, NX, NY, NXY, and oxide data 

-X and Y: x-y coordinates EPMA position, NX and NY: pixel number in x-y (matrix)
-NXY: consecutive pixel number, all the rest values are weight percents including the total.
-Oxide data should be reported in WT% (e.g. SiO2 WT%, TiO2 WT%, Al2O3 WT%, MgO WT%

Regarding single-points analyses, the heading of the file should contain at least sample name and the oxides"""

show_intro(Introd)

###########################################################################################################################################################
# Loading EPMA data
###########################################################################################################################################################
attempts = 3  #Number of attempts allowed to input file's path

for attempt in range(attempts):
    file_path = input('Please type the path where the data file is located:   ') #Message to get input from user
    #changing name of the file
    extraname = r'_FirstPass_\d{5}__Oxide_Image_Classify' #Part of the file name
    new_file_name = re.sub(extraname, '', os.path.basename(file_path)) #simplifies the file name by removing the previous part
    directory = os.path.dirname(file_path) #directory of the input excluding the name
    new_file_path = os.path.join(directory, new_file_name) #new path using the new name
    os.rename(file_path, new_file_path) #renames the original file
    print(f"File renamed to: {new_file_path}")

    with open (new_file_path, 'r') as file:
        lines = file.readlines()
    lines[1] = lines[1].replace('\t\t', '\t')

    with open(new_file_path, 'w') as corrected_file:
        corrected_file.writelines(lines)

    startcounting = time.time() #timer starts
    try:
        data = pd.read_csv(new_file_path, sep="\t", header=1)  #Reading data file from specified path but skipping the first row of data
        print(f'Loading of data was successful! Total elapsed time: {time.time() - startcounting:.1f} seconds')
        break  #exit the loop if file exists
    except FileNotFoundError:
        #goes here if file cannot be read or does not exist
        print(f"File not found. Please check the path and try again. Attempt {attempt + 1} of {attempts}.")
        if attempt == attempts - 1:
            #Exits the software after two attempts
            print("Failed to load the file after multiple attempts. Exiting epMapy.")
            raise SystemExit
print(data.head(5))
###########################################################################################################################################################
# Deciding how to proceed
###########################################################################################################################################################
#map_point = input('Do you want to process a map: (Yes/No)   ') #Message to get input from user
#if map_point == 'yes':
        #Asks user which oxide to use for selecting the profile points
        #selected_oxide = input(f"Which oxide would you like to use for selecting the profile points? Choose from: {', '.join(selected_oxides)}: ").strip()

        #if selected_oxide not in selected_oxides:
            #print("Invalid oxide selection. Please select an oxide from the list.")
            #continue
        #extract_profile(oxide_grids,selected_oxides, pixel_size, selected_oxide, sample_name)


###########################################################################################################################################################
# Initialising and filling dictionaries for compositions and coordinates.
###########################################################################################################################################################
oxide_columns = {} 
coord_columns = {}
oxide_columns_original = {}

#Removes 'WT%' from column names but keeps capitalisation
data.columns = data.columns.str.replace(" WT%", "", regex=True).str.replace("wt%", "", regex=True).str.strip()

#Iterate over column names for matching with potential oxides and coordinates
for norm_col in data.columns:
    if re.search(r"(SiO2|TiO2|Cr2O3|Al2O3|FeO|MnO|NiO|MgO|CaO|Na2O|K2O|P2O5|F|Cl|SO3|BaO|SrO|La2O3|Ce2O3|Total)", norm_col, re.IGNORECASE):
        print(f"Matched oxide column: {norm_col}") #prints which oxides are in the file
        oxide_columns[norm_col] = data[norm_col].values  #Fills in the values of each oxide in the dictionary
        oxide_columns_original[norm_col] = data[norm_col].values  #clone dictionary to keep original column names for following steps
    elif re.search(r"(X|Y|NX|NY)", norm_col, re.IGNORECASE):
        print(f"Matched coordinate column: {norm_col}")
        coord_columns[norm_col] = data[norm_col].values  #Fills in the values of each coordinate or pixel number in the dictionary

#List of potential oxides measured with the microprobe
expected_oxides = ["SiO2", "TiO2", "Cr2O3", "Al2O3", "FeO", "MnO", "NiO", "MgO", "CaO", 
                   "Na2O", "K2O", "P2O5", "F", "Cl", "SO3", "BaO", "SrO", "La2O3", "Ce2O3", "Total"]
#Fills missing oxide columns with very small numbers to avoid dividing by zero or NAN.
for oxide in expected_oxides:
    if oxide not in oxide_columns:
        print(f"{oxide} column not found. Filling with 0.00000000001 to avoid issues when doing calculations.")
        oxide_columns[oxide] = np.full(len(data), 0.00000000001)  #Fills with 0.00000000001

###########################################################################################################################################################
# Extract column heads and convert dictionaries into numpy arrays
###########################################################################################################################################################
#The following dictionary also includes also those oxides that were filled with 0.00000000001 for calculation purposes
oxide_column_indices = {col: index for index, col in enumerate(oxide_columns.keys())}
#The following clone dictiornary doesnt contain the oxides that were filled with 0.000000000001 values
oxide_column_indices_original = {col: index for index, col in enumerate(oxide_columns_original.keys())}
#coordinates and pixel number dictiornary
coord_column_indices = {col: index for index, col in enumerate(coord_columns.keys())}
#numpy arrays
oxide_data = np.array([oxide_columns[key] for key in oxide_columns])
coord_data = np.array([coord_columns[key] for key in coord_columns])

###########################################################################################################################################################
#Structuring data so that it can be plotted
###########################################################################################################################################################
c, d = int(max(coord_data[coord_column_indices["NX"]])), int(max(coord_data[coord_column_indices["NY"]])) #Grid size matches the maximum NX and NY values

# Initialize arrays for various properties
# molar indexes: NK_A:Na2O+K2O/Al2O3, NKC_A:Na2O+K2O+CaO/Al2O3, K_Na:K2O/Na2O, M_f:
NK_A, NKC_A, Mg_MgFe2, K_Na, M_f = [np.zeros((c, d)) for _ in range(5)]

#creating grids with coordinates and
coord_grids = {coordinate: np.zeros((c,d)) for coordinate in coord_column_indices.keys()}
oxide_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_original.keys()}
oxidef_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_original.keys()}
data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
print(oxide_column_indices)
fe2o3_included, oxide_column_indices_fe = add_fe2o3(data_majors,oxide_column_indices) #calculates fe2o3 and updates the oxide_column_indices
print(oxide_column_indices_fe,oxide_column_indices)
oxidef_grids_fe = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices_fe.keys()}
data_anhf = to_anhydrous(data_majors, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
data_anhf_fe = to_anhydrous(fe2o3_included, oxide_column_indices_fe) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
data_molf_fe   = to_mol(data_anhf_fe, oxide_column_indices_fe)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
data_mol   = to_mol(data_majors, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
data_cat   = to_cat(data_mol, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
data_normf  = norm_calc(data_molf_fe, oxide_column_indices_fe) #calculates normative mineralogy

#oxidefe3_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}

#molanh_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}

# Loop through all data points to fill grid arrays
for idx in range(len(coord_data[coord_column_indices["NXY"]])):
    nx, ny = int(coord_data[coord_column_indices["NX"]][idx]) - 1, int(coord_data[coord_column_indices["NY"]][idx]) - 1
    
    if data_majors[oxide_column_indices["Total"]][idx] > 0:
        for oxide, col_index in oxide_column_indices_original.items():
            col_index = oxide_column_indices[oxide]
            oxide_grids[oxide][nx, ny] = data_majors[col_index][idx]
            oxidef_grids[oxide][nx, ny] = data_anhf[col_index][idx]
            oxidef_grids_fe[oxide][nx, ny] = data_anhf_fe[col_index][idx]
        for coordinate, col_index in coord_column_indices.items():
            col_index = coord_column_indices[coordinate]
            coord_grids[coordinate][nx, ny] = coord_data[col_index][idx]

        # Fill derived properties
        M_f[nx, ny] = (
            (data_catf[oxide_column_indices["Na2O"],idx] +
             data_catf[oxide_column_indices["K2O"],idx] +
             (data_catf[oxide_column_indices["CaO"],idx] * 2)) /
            (data_catf[oxide_column_indices["SiO2"],idx] *
             data_catf[oxide_column_indices["Al2O3"],idx])
        )
        NK_A[nx, ny] = (
            (data_molf[oxide_column_indices["Na2O"],idx] +
             data_molf[oxide_column_indices["K2O"],idx]) /
            data_molf[oxide_column_indices["Al2O3"],idx]
        )
        NKC_A[nx, ny] = (
            (data_molf[oxide_column_indices["Na2O"],idx] +
             data_molf[oxide_column_indices["K2O"],idx] +
             data_molf[oxide_column_indices["CaO"],idx]) /
            data_molf[oxide_column_indices["Al2O3"],idx]
        )
        K_Na[nx, ny] = (
            data_anhf[oxide_column_indices["K2O"],idx] * 8302 /
            data_anhf[oxide_column_indices["Na2O"],idx] * 7419
        )
    else:
        for oxide in oxide_column_indices_original.keys():
            oxide_grids[oxide][nx, ny] = -1  # Placeholder for invalid points
            oxidef_grids[oxide][nx, ny] = -1  # Placeholder for invalid points
        for grid in [M_f, NK_A, NKC_A, K_Na]:
            grid[nx, ny] = -1
        for coordinate, col_index in coord_column_indices.items():
            col_index = coord_column_indices[coordinate]
            coord_grids[coordinate][nx, ny] = coord_data[col_index][idx]


###########################################################################################################################################################
# Plotting and profile extraction loop
###########################################################################################################################################################
# Main loop for plotting and profile extraction
sample_name = input("Enter the sample name: ").strip()
pixel_size = (coord_data[coord_column_indices["X"], 1] - coord_data[coord_column_indices["X"], 0]) * 1000

while True:
    selected_oxides_input = input("Enter the oxides you want to plot, separated by commas, or type 'all' to plot all: ").strip()
    if selected_oxides_input.lower() == 'all':
        selected_oxides = list(oxide_columns_original.keys())  # Select all available oxides
    else:
        selected_oxides = [oxide.strip() for oxide in selected_oxides_input.split(',')]  # Creates a list with the selected oxides

    cm = 1 / 2.54
    plt.figure(figsize=(6.5 * cm, 6.5 * cm))
    fig, axs = plt.subplots(1, len(selected_oxides), figsize=(3 * len(selected_oxides), 3))

    for i, oxide in enumerate(selected_oxides):
        oxide = oxide.strip()
        if oxide in oxide_grids:
            ax = axs[i] if len(selected_oxides) > 1 else axs
            im = ax.imshow(oxide_grids[oxide], cmap='viridis', vmin=0, vmax=np.max(oxide_grids[oxide]))
            im.cmap.set_under('black')
            ax.set_title(f"{oxide} wt.%")
            ax.axis("off")
            fig.colorbar(im, ax=ax)
        else:
            print(f"Warning: {oxide} is not a valid oxide.")

    scalebar = ScaleBar(pixel_size, "um")
    plt.gca().add_artist(scalebar)
    plt.tight_layout()

    formatted_oxides = "_".join([oxide.strip() for oxide in selected_oxides])
    pdf_filename = f"{sample_name}_{formatted_oxides}.pdf"
    plt.savefig(pdf_filename, dpi=600, transparent=True, bbox_inches='tight')
    print(f"Plot saved as {pdf_filename}")
    plt.show()

    # Now, enter a loop to allow multiple profile extractions
    profile_index = 1  # Initialize profile index
    while True:
        # Ask User if they want to extract a profile
        extract = input("Would you like to extract a traverse for the selected oxides? (yes/no): ").strip().lower()
        if extract == 'yes':
            # Ask user which oxide to use for selecting the profile points
            selected_oxide = input(f"Which oxide would you like to use for selecting the profile points? Choose from: {', '.join(selected_oxides)}: ").strip()

            if selected_oxide not in selected_oxides:
                print("Invalid oxide selection. Please select an oxide from the list.")
                continue

            # Extract profile for the selected oxide
            extract_profile(oxide_grids, selected_oxides, pixel_size, selected_oxide, sample_name, profile_index)
            profile_index += 1  # Increment the profile index after each extraction

        else:
            break  # Exit the profile extraction loop if the user is done

    # Ask if the user wants to plot another set of oxides after extracting profiles
    repeat = input("Would you like to plot another set of oxides? (yes/no): ").strip().lower()
    if repeat != 'yes':
        break  # Exit the entire loop if the user doesn't want to plot another set
