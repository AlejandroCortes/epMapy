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
d, c = int((max(coord_data[:,2]))-1), int((max(coord_data[:,3]))-1) # Determine grid size from unique coordinates
NK_A,NKC_A,Mg_MgFe2,K_Na                                        = [np.zeros((d,c)) for i in range(4)] #setting up empty molar fractions arrays
Mf                                                              = np.zeros((d,c)) #setting up empty cation fractions arrays (M factor)
SiO2p,TiO2p,Al2O3p,FeOp,MnOp,MgOp,CaOp,Na2Op,K2Op,P2O5p,totp    = [np.zeros((d,c)) for i in range(11)] #setting up empty oxides wt.% arrays
SiO2,TiO2,Al2O3,FeO,MnO,MgO,CaO,Na2O,K2O,P2O5,tot               = [np.zeros((d,c)) for i in range(11)] #setting up empty anhydrous-based oxides wt.% arrays
DI                                                              = np.zeros((d,c)) #setting up empty Normative Differentiation Index array
Q,ort,lc,ab,ne,an,C,ac,ns,di,ol,hy,mt,ilm,ap,Fsp,mrg,px,ox      = [np.zeros((d,c)) for i in range(19)] #setting up empty normative minerals wt.% arrays
Nx,Ny                                                           = [np.zeros((d,c)) for i in range(2)] #setting up empty coordinate arrays (mesh-like)

#filtered_df = df[df['emails'].str.contains('GMAIL', case=False)]
    #reorganising calculated data based on pixel coordinates
    #filling up empty arrays
for j in range (d):
    for i in range (c):
        if data_majors[j*(c+1)+i,-1]>0:
            SiO2p[j,i],TiO2p[j,i],Al2O3p[j,i],FeOp[j,i],MnOp[j,i],MgOp[j,i],CaOp[j,i],Na2Op[j,i],K2Op[j,i],P2O5p[j,i],totp[j,i]     = [data_majors[j*c+j+i,q] for q in range (11)]
            SiO2[j,i],TiO2[j,i],Al2O3[j,i],FeO[j,i],MnO[j,i],MgO[j,i],CaO[j,i],Na2O[j,i],K2O[j,i],P2O5[j,i],tot[j,i]                = [data_anhf[j*c+j+i,q] for q in range (11)]
            #Q[j,i],ort[j,i],lc[j,i],ab[j,i],ne[j,i],an[j,i],C[j,i],ac[j,i],ns[j,i],di[j,i],ol[j,i],hy[j,i],mt[j,i],ilm[j,i],ap[j,i] = [data_normf[j*c+j+i,q] for q in range (15)]
            Mf[j,i]        = (data_catf[j*(c+1)+i,7]+data_catf[j*(c+1)+i,8]+
                            (data_catf[j*(c+1)+i,6]*2))/(data_catf[j*(c+1)+i,0]*data_catf[j*(c+1)+i,2])
            NK_A[j,i]      = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8])/data_molf[j*(c+1)+i,2]
            NKC_A[j,i]     = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8]+
                             data_molf[j*(c+1)+i,6])/data_molf[j*(c+1)+i,2]
            #Mg_MgFe2[j,i]  = data_molf2[j*(c+1)+i,6]/(data_molf2[j*(c+1)+i,4]+data_molf2[j*(c+1)+i,6]+
            #                data_molf2[j*(c+1)+i,5])
            K_Na[j,i]      = (data_anhf[j*c+j+i,8]*8302)/(data_anhf[j*c+j+i,7]*7419)
            #DI[j,i]        = np.sum(data_normf[j*c+j+i,0:5])+data_normf[j*c+j+i,8]
            #Fsp[j,i]       = ort[j,i]+ab[j,i]+an[j,i]
            #px[j,i]        = ac[j,i]+di[j,i]+hy[j,i]
            #ox[j,i]        = mt[j,i]+ilm[j,i]
            #Nx[j,i]        = data_coord[j*c+j+i,2]
            #Ny[j,i]        = data_coord[j*c+j+i,3]
        else:
            SiO2p[j,i],TiO2p[j,i],Al2O3p[j,i],FeOp[j,i],MnOp[j,i],MgOp[j,i],CaOp[j,i],Na2Op[j,i],K2Op[j,i],P2O5p[j,i]               = [300 for q in range (10)]
            SiO2[j,i],TiO2[j,i],Al2O3[j,i],FeO[j,i],MnO[j,i],MgO[j,i],CaO[j,i],Na2O[j,i],K2O[j,i],P2O5[j,i]                         = [300 for q in range (10)]
            totp[j,i],tot[j,i],Mf[j,i],K_Na[j,i],Q[j,i]                                                                             = [-1 for q in range (5)]
            NK_A[j,i],NKC_A[j,i],Mg_MgFe2[j,i],DI[j,i],Fsp[j,i],px[j,i],ox[j,i]                                                     = [300 for q in range (7)]
            #ort[j,i],lc[j,i],ab[j,i],ne[j,i],an[j,i],C[j,i],ac[j,i],ns[j,i],di[j,i],ol[j,i],hy[j,i],mt[j,i],ilm[j,i],ap[j,i]        = [300 for q in range (14)]
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
plt.savefig(filename+'_oxides_anh.pdf',dpi=600, transparent=True, bbox_inches='tight')