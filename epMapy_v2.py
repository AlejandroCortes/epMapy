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
from functionsmapping import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data  # importing external functions from functionsmapping.py

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

for col in data.columns: # Identify oxide and coordinate columns by matching strings
    if re.search(r"(SiO2|TiO2|Cr2O3|Al2O3|FeO|MnO|NiO|MgO|CaO|Na2O|K2O|P2O5|F|Cl|SO3|BaO|SrO|Total)", col, re.IGNORECASE): # Identify oxides by molecular/element formula
        oxide_columns[col] = data[col].values  # filling the keys with available oxides
    elif re.search(r"(X|Y|NX|NY)", col, re.IGNORECASE) and "WT%" not in col:# Identify spatial or pyxel coordinates by string match
        coord_columns[col] = data[col].values # filling the keys: 'X','Y','NX','NY'

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

###########################################################################################################################################################
# Calling multiple functions
###########################################################################################################################################################
data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
data_anhf = to_anhydrous(oxide_data, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf)         #calculates anhydrous-based cation compositions mol
data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors))) #calculates oxide compositions including Fe2O3 mol
data_normf  = norm_calc(data_molf2) #calculates normative mineralogy

###########################################################################################################################################################
# Structuring data so that it can be plotted
###########################################################################################################################################################
d, c = int((max(coord_data[:,2]))-1), int((max(coord_data[:,3]))-1) # Determine grid size from unique coordinates
NK_A,NKC_A,Mg_MgFe2,K_Na                                        = [np.zeros((d,c)) for i in range(4)] #setting up empty molar fractions arrays
Mf                                                              = np.zeros((d,c)) #setting up empty cation fractions arrays (M factor)
SiO2p,TiO2p,Al2O3p,FeOp,MnOp,MgOp,CaOp,Na2Op,K2Op,P2O5p,totp    = [np.zeros((d,c)) for i in range(11)] #setting up empty oxides wt.% arrays
SiO2,TiO2,Al2O3,FeO,MnO,MgO,CaO,Na2O,K2O,P2O5,tot               = [np.zeros((d,c)) for i in range(11)] #setting up empty anhydrous-based oxides wt.% arrays
DI                                                              = np.zeros((d,c)) #setting up empty Normative Differentiation Index array
Q,ort,lc,ab,ne,an,C,ac,ns,di,ol,hy,mt,ilm,ap,Fsp,mrg,px,ox      = [np.zeros((d,c)) for i in range(19)] #setting up empty normative minerals wt.% arrays
Nx,Ny                                                           = [np.zeros((d,c)) for i in range(2)] #setting up empty coordinate arrays (mesh-like)
