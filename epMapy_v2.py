# Import libraries
import numpy as np
import pandas as pd
import re  # for pattern matching
import time

from functionsmapping import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data #importing external functions

# Load data
name = input('Please type the path where the data file is located:   ')
try:
    data = pd.read_excel(name)
    print(f'Loading of data was successful! Total elapsed time: {time.time() - st0:.1f} seconds')
except FileNotFoundError:
    print(f"File not found. Total execution time: {time.time() - st0:.1f} seconds")
    raise SystemExit

# Detect columns based on patterns
oxide_columns = {}
coord_columns = {}

# Define expected oxides for later filling missing ones
expected_oxides = ["SiO2", "TiO2", "Cr2O3", "Al2O3", "FeO", "MnO", "NiO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "F", "Cl", "SO3", "BaO", "SrO", "Total"]

# Identify oxide and coordinate columns based on patterns
for col in data.columns:
    # Identify oxides by keywords
    if re.search(r"(SiO2|TiO2|Cr2O3|Al2O3|FeO|MnO|NiO|MgO|CaO|Na2O|K2O|P2O5|F|Cl|SO3|BaO|SrO|Total)", col, re.IGNORECASE):
        oxide_columns[col] = data[col].fillna(0.000000000000000001).values  # Fill NaN with 0 for consistency
    elif re.search(r"(X|Y|NX|NY)", col, re.IGNORECASE) and "WT%" not in col:
        coord_columns[col] = data[col].values

# Fill missing oxide columns with small numbers
for oxide in expected_oxides:
    if oxide not in oxide_columns:
        print(f"{oxide} column not found. Filling with 0.00000000001.")
        oxide_columns[oxide] = np.full(len(data), 0.00000000001)  # Fill with 0.00000000001

# Convert oxide and coordinate data to NumPy arrays
oxide_column_indices = {col: idx for idx, col in enumerate(oxide_columns.keys())}
oxide_data = np.array([oxide_columns[key] for key in oxide_columns])
coord_data = np.array([coord_columns[key] for key in coord_columns])

# Determine grid size from unique coordinates
d, c = int((max(coord_data[:,2]))-1), int((max(coord_data[:,3]))-1)

#setting up empty molar fractions arrays
NK_A,NKC_A,Mg_MgFe2,K_Na                                        = [np.zeros((d,c)) for i in range(4)]
#setting up empty cation fractions arrays
Mf                                                              = np.zeros((d,c)) #M factor
#setting up empty oxides wt.% arrays
SiO2p,TiO2p,Al2O3p,FeOp,MnOp,MgOp,CaOp,Na2Op,K2Op,P2O5p,totp    = [np.zeros((d,c)) for i in range(11)]
#setting up empty anhydrous-based oxides wt.% arrays
SiO2,TiO2,Al2O3,FeO,MnO,MgO,CaO,Na2O,K2O,P2O5,tot               = [np.zeros((d,c)) for i in range(11)]
#setting up empty Normative Differentiation Index array
DI                                                              = np.zeros((d,c))
#setting up empty normative minerals wt.% arrays
Q,ort,lc,ab,ne,an,C,ac,ns,di,ol,hy,mt,ilm,ap,Fsp,mrg,px,ox      = [np.zeros((d,c)) for i in range(19)]
#setting up empty coordinate arrays (mesh-like)
Nx,Ny                                                           = [np.zeros((d,c)) for i in range(2)] #coordinates of pixels (columns)


#employing functions
data_majors = clean_data(oxide_data, oxide_column_indices)   #cleaning holes and cracks that retrive very low totals
data_majors = clean_data(oxide_data, oxide_column_indices) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf)         #calculates anhydrous-based cation compositions mol
data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors))) #calculates oxide compositions including Fe2O3 mol
data_normf  = norm_calc(data_molf2) #calculates normative mineralogy
