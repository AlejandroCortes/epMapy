introd = f'''
This script allows plotting EPMA quantitative maps and performing calculations with the data within, as well as treatment of single-point data acquired with the probe.
Please be aware that this script is designed for the output file obtained with CalcImage (J. Donovan).
The input file for this script is a .xlsx file, so you should convert .DAT to .xlsx before
using it. It must contain at least the following column names (case and character sensitive):
X, Y, NX, NY, NXY, SiO2 WT%, TiO2 WT%, Al2O3 WT%, FeO WT%, MnO WT%, MgO WT%, CaO WT%, Na2O WT%, K2O WT%, P2O5 WT%, Total
X and Y: x-y coordinates EPMA position, NX and NY: pixel number in x-y (matrix)
NXY: consecutive pixel number, all the rest values are weight percents including the total.
Unfortunately, this script does not work with more oxides or halogens at the moment, but can be easily
modified to do so by the user. Just keep in mind to also modify the external functions that are being used.
'''
print(introd)

###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
import numpy as np  # to deal mathematically with multidimensional array objects
import pandas as pd  # to structure data sets
import re  # to compare a set of string-type attributes
import time  # to handle time-related calculations
import matplotlib.pyplot as plt #default plotting library
from functionsmapping import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data  # importing external functions from functionsmapping.py
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
import tkinter
###########################################################################################################################################################
# Loading EPMA data
###########################################################################################################################################################
attempts = 2  # number of attempts allowed to input file's path

for attempt in range(attempts):
    file_path = input('Please type the path where the data file is located:   ')
    st0 = time.time()
    try:
        data = pd.read_excel(file_path)  # try reading data file from specified path
        print(f'Loading of data was successful! Total elapsed time: {time.time() - st0:.1f} seconds')
        break  # if successful, exit the loop
    except FileNotFoundError:
        print(f"File not found. Please check the path and try again. Attempt {attempt + 1} of {attempts}.")
        if attempt == attempts - 1:
            print("Failed to load the file after multiple attempts. Exiting the software.")
            raise SystemExit

###########################################################################################################################################################
# Initialising and filling dictionaries for compositions and coordinates.
###########################################################################################################################################################
oxide_columns = {} 
coord_columns = {}
oxide_columns_original = {}

# Normalize column names by removing 'WT%' but retain capitalization
data.columns = data.columns.str.replace(" WT%", "", regex=True).str.replace("wt%", "", regex=True).str.strip()
# Create a mapping of original column names to normalized column names
#column_map = dict(zip(normalized_columns, data.columns))

# Iterate over normalized column names for matching
for norm_col in data.columns:
    if re.search(r"(SiO2|TiO2|Cr2O3|Al2O3|FeO|MnO|NiO|MgO|CaO|Na2O|K2O|P2O5|F|Cl|SO3|BaO|SrO|Total)", norm_col, re.IGNORECASE):
        print(f"Matched oxide column: {norm_col}")
        oxide_columns[norm_col] = data[norm_col].values  # Use the original column name
        oxide_columns_original[norm_col] = data[norm_col].values  # Use the original column name
    elif re.search(r"(X|Y|NX|NY)", norm_col, re.IGNORECASE):
        print(f"Matched coordinate column: {norm_col}")
        coord_columns[norm_col] = data[norm_col].values  # Use the original column name

# Define potential oxides or elements that will be needed for calculations
expected_oxides = ["SiO2", "TiO2", "Cr2O3", "Al2O3", "FeO", "MnO", "NiO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "F", "Cl", "SO3", "BaO", "SrO", "Total"]

# Fill missing oxide columns with very small numbers
for oxide in expected_oxides:
    if oxide not in oxide_columns:
        print(f"{oxide} column not found. Filling with 0.00000000001 to avoid issues when doing calculations.")
        oxide_columns[oxide] = np.full(len(data), 0.00000000001)  # Fill with 0.00000000001

###########################################################################################################################################################
# Convert dictionaries into numpy arrays but keeping their column indexes
###########################################################################################################################################################
oxide_column_indices = {col: index for index, col in enumerate(oxide_columns.keys())}
oxide_column_indices_original = {col: index for index, col in enumerate(oxide_columns_original.keys())}
coord_column_indices = {col: index for index, col in enumerate(coord_columns.keys())}
oxide_data = np.array([oxide_columns[key] for key in oxide_columns])
coord_data = np.array([coord_columns[key] for key in coord_columns])


###########################################################################################################################################################
# Calling multiple functions
###########################################################################################################################################################
data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
data_anhf = to_anhydrous(oxide_data, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
fe2o3_included, oxide_column_indices = add_fe2o3(data_majors,oxide_column_indices)
fe2o3_anh = to_anhydrous(fe2o3_included, oxide_column_indices)
fe2o3_mol = to_mol(fe2o3_anh, oxide_column_indices)
#data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors, oxide_column_indices), oxide_column_indices),oxide_column_indices) #calculates oxide compositions including Fe2O3 mol
data_normf  = norm_calc(fe2o3_mol, oxide_column_indices) #calculates normative mineralogy

###########################################################################################################################################################
# Structuring data so that it can be plotted
###########################################################################################################################################################
# Ensure the grid size matches the maximum NX and NY values
c, d = int(max(coord_data[coord_column_indices["NX"]])), int(max(coord_data[coord_column_indices["NY"]]))

# Initialize arrays for various properties
NK_A, NKC_A, Mg_MgFe2, K_Na = [np.zeros((c, d)) for _ in range(4)]
Mf = np.zeros((c, d))
oxide_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}
oxideanh_grids = {oxide: np.zeros((c, d)) for oxide in oxide_column_indices.keys()}

# Loop through all data points to fill grid arrays
for idx in range(len(coord_data[coord_column_indices["NXY"]])):
    nx, ny = int(coord_data[coord_column_indices["NX"]][idx]) - 1, int(coord_data[coord_column_indices["NY"]][idx]) - 1
    
    if data_majors[oxide_column_indices["Total"]][idx] > 0:
        for oxide, col_index in oxide_column_indices_original.items():
            col_index = oxide_column_indices[oxide]
            oxide_grids[oxide][nx, ny] = data_majors[col_index][idx]
            oxideanh_grids[oxide][nx, ny] = data_anhf[col_index][idx]

        # Fill derived properties
        Mf[nx, ny] = (
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
            oxide_grids[oxide][nx, ny] = 300  # Placeholder for invalid points
            oxideanh_grids[oxide][nx, ny] = 300  # Placeholder for invalid points
        for grid in [Mf, NK_A, NKC_A, K_Na]:
            grid[nx, ny] = -1


# Prompt user to select oxides to plot
print("Available oxides to plot:", list(oxide_columns_original.keys()))
selected_oxides = input("Enter the oxides you want to plot, separated by commas: ").split(',')

# Plot the selected oxides
cm = 1 / 2.54
plt.figure(figsize=(6.5 * cm, 6.5 * cm))
fig, axs = plt.subplots(1, len(selected_oxides), figsize=(3 * len(selected_oxides), 3))

for i, oxide in enumerate(selected_oxides):
    oxide = oxide.strip()  # Remove extra spaces
    if oxide in oxide_grids:
        ax = axs[i] if len(selected_oxides) > 1 else axs
        im = ax.imshow(oxide_grids[oxide], cmap='viridis', vmin=0, vmax=100)
        im.cmap.set_under('black')
        im.cmap.set_over('black')
        ax.set_title(f"{oxide} wt.%")
        ax.axis("off")
        fig.colorbar(im, ax=ax)
    else:
        print(f"Warning: {oxide} is not a valid oxide.")

scalebar = ScaleBar(0.000005)  # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig('selected_oxides.pdf', dpi=600, transparent=True, bbox_inches='tight')
plt.show()