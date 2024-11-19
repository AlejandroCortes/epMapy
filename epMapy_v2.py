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

# Normalize column names by removing 'WT%' but retain capitalization
data.columns = data.columns.str.replace(" WT%", "", regex=True).str.replace("wt%", "", regex=True).str.strip()
# Create a mapping of original column names to normalized column names
#column_map = dict(zip(normalized_columns, data.columns))

# Iterate over normalized column names for matching
for norm_col in data.columns:
    if re.search(r"(SiO2|TiO2|Cr2O3|Al2O3|FeO|MnO|NiO|MgO|CaO|Na2O|K2O|P2O5|F|Cl|SO3|BaO|SrO|Total)", norm_col, re.IGNORECASE):
        print(f"Matched oxide column: {norm_col}")
        oxide_columns[norm_col] = data[norm_col].values  # Use the original column name
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
oxide_data = np.array([oxide_columns[key] for key in oxide_columns])
coord_data = np.array([coord_columns[key] for key in coord_columns])
print(oxide_data[oxide_column_indices["SiO2"]])
###########################################################################################################################################################
# Calling multiple functions
###########################################################################################################################################################
data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
data_anhf = to_anhydrous(oxide_data, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf, oxide_column_indices)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf, oxide_column_indices)         #calculates anhydrous-based cation compositions mol
#data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors, oxide_column_indices), oxide_column_indices),oxide_column_indices) #calculates oxide compositions including Fe2O3 mol
#data_normf  = norm_calc(data_molf2, oxide_column_indices) #calculates normative mineralogy

###########################################################################################################################################################
# Structuring data so that it can be plotted
###########################################################################################################################################################
# Ensure the grid size matches the maximum NX and NY values
d, c = int(max(coord_columns["NX"])), int(max(coord_columns["NY"]))

# Initialize arrays for various properties
NK_A, NKC_A, Mg_MgFe2, K_Na = [np.zeros((d, c)) for _ in range(4)]
Mf = np.zeros((d, c))
SiO2p, TiO2p, Al2O3p, FeOp, MnOp, MgOp, CaOp, Na2Op, K2Op, P2O5p, totp = [np.zeros((d, c)) for _ in range(11)]
SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5, tot = [np.zeros((d, c)) for _ in range(11)]

# Loop through all data points to fill grid arrays
for idx in range(len(data_majors)):
    # Convert pixel number to grid indices
    nx, ny = int(coord_columns["NX"][idx]) - 1, int(coord_columns["NY"][idx]) - 1

    # Check if Total column value is valid
    if data_majors[idx, oxide_column_indices["Total"]] > 0:
        # Fill oxide weight percent grids using column names
        SiO2p[nx, ny] = data_majors[idx, oxide_column_indices["SiO2"]]
        TiO2p[nx, ny] = data_majors[idx, oxide_column_indices["TiO2"]]
        Al2O3p[nx, ny] = data_majors[idx, oxide_column_indices["Al2O3"]]
        FeOp[nx, ny] = data_majors[idx, oxide_column_indices["FeO"]]
        MnOp[nx, ny] = data_majors[idx, oxide_column_indices["MnO"]]
        MgOp[nx, ny] = data_majors[idx, oxide_column_indices["MgO"]]
        CaOp[nx, ny] = data_majors[idx, oxide_column_indices["CaO"]]
        Na2Op[nx, ny] = data_majors[idx, oxide_column_indices["Na2O"]]
        K2Op[nx, ny] = data_majors[idx, oxide_column_indices["K2O"]]
        P2O5p[nx, ny] = data_majors[idx, oxide_column_indices["P2O5"]]
        totp[nx, ny] = data_majors[idx, oxide_column_indices["Total"]]

        # Fill anhydrous oxide grids
        SiO2[nx, ny] = data_anhf[idx, oxide_column_indices["SiO2"]]
        TiO2[nx, ny] = data_anhf[idx, oxide_column_indices["TiO2"]]
        Al2O3[nx, ny] = data_anhf[idx, oxide_column_indices["Al2O3"]]
        FeO[nx, ny] = data_anhf[idx, oxide_column_indices["FeO"]]
        MnO[nx, ny] = data_anhf[idx, oxide_column_indices["MnO"]]
        MgO[nx, ny] = data_anhf[idx, oxide_column_indices["MgO"]]
        CaO[nx, ny] = data_anhf[idx, oxide_column_indices["CaO"]]
        Na2O[nx, ny] = data_anhf[idx, oxide_column_indices["Na2O"]]
        K2O[nx, ny] = data_anhf[idx, oxide_column_indices["K2O"]]
        P2O5[nx, ny] = data_anhf[idx, oxide_column_indices["P2O5"]]
        tot[nx, ny] = data_anhf[idx, oxide_column_indices["Total"]]

        # Fill derived properties
        Mf[nx, ny] = (
            (data_catf[idx, oxide_column_indices["Na2O"]] +
             data_catf[idx, oxide_column_indices["K2O"]] +
             (data_catf[idx, oxide_column_indices["CaO"]] * 2)) /
            (data_catf[idx, oxide_column_indices["SiO2"]] *
             data_catf[idx, oxide_column_indices["Al2O3"]])
        )
        NK_A[nx, ny] = (
            (data_molf[idx, oxide_column_indices["Na2O"]] +
             data_molf[idx, oxide_column_indices["K2O"]]) /
            data_molf[idx, oxide_column_indices["Al2O3"]]
        )
        NKC_A[nx, ny] = (
            (data_molf[idx, oxide_column_indices["Na2O"]] +
             data_molf[idx, oxide_column_indices["K2O"]] +
             data_molf[idx, oxide_column_indices["CaO"]]) /
            data_molf[idx, oxide_column_indices["Al2O3"]]
        )
        K_Na[nx, ny] = (
            data_anhf[idx, oxide_column_indices["K2O"]] * 8302 /
            data_anhf[idx, oxide_column_indices["Na2O"]] * 7419
        )
    else:
        # Fill invalid points with placeholder values
        for grid in [SiO2p, TiO2p, Al2O3p, FeOp, MnOp, MgOp, CaOp, Na2Op, K2Op, P2O5p, SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5]:
            grid[nx, ny] = 300
        for grid in [totp, tot, Mf, NK_A, NKC_A, K_Na]:
            grid[nx, ny] = -1


cm       = 1/2.54

plt.figure(figsize=(6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(SiO2,cmap='viridis',vmin=0, vmax=100)
a0.cmap.set_over('black')
axs[0,0].set_title("SiO$_{2}$ wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Al2O3,cmap='viridis',vmin=0, vmax=45)
a1.cmap.set_over('black')
axs[0,1].set_title("Al$_{2}$O$_{3}$ wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(FeO,cmap='viridis',vmin=0, vmax=100)
a2.cmap.set_over('black')
axs[0,2].set_title("FeO wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(MgO,cmap='viridis',vmin=0, vmax=45)
a3.cmap.set_over('black')
axs[1,0].set_title("MgO wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(Na2O,cmap='viridis',vmin=0, vmax=15)
a4.cmap.set_over('black')
axs[1,1].set_title("Na$_{2}$O wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(P2O5,cmap='viridis',vmin=0, vmax=60)
a5.cmap.set_over('black')
axs[1,2].set_title("P$_{2}$O$_{5}$")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig('oxides_anh.pdf',dpi=600, transparent=True, bbox_inches='tight')