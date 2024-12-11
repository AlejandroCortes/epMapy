#Importing libraries a$nd changing plotting parameters

import numpy as np #math functions as pi,e,random,arrays, etc 
import matplotlib.pyplot as plt #default plotting library
import pandas as pd #for printing nice tables (pd.DataFrame)
import seaborn as sbn #extra version of plt
plt.rcParams['pdf.fonttype'] = 42 #to use Type 42 (a.k.a. TrueType) fonts for PDF files
plt.rcParams['ps.fonttype']  = 42 #to use Type 42 (a.k.a. TrueType) fonts for PostScript files
plt.rcParams['svg.fonttype'] = 'none' #for neither embedding the font nor rendering the text as path
plt.rcParams.update({'font.size': 7}) #defined font size
plt.rcParams['figure.dpi']   = 200 #increases the quality of images
import mpltern #to plot ternary diagrams
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
from functionsepmapy import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3 #importing external functions

print('\n\nThis script allows plotting EPMA quantitative maps and do calculations with the data within.\nPlease be aware that this script is designed for the output file obtained with CalcImage (J. Donovan).\nThe input file for this script is a .xslx file, so you should convert .DAT to .xslx before\nat least the following column names (case and character sensitive):\nX,Y,NX,NY,NXY,SiO2 WT%,TiO2 WT%,Al2O3 WT%,FeO WT%,MnO WT%,MgO WT%,CaO WT%,Na2O WT%,K2O WT%,P2O5 WT%,Total\nX and Y: x-y coordinates EPMA position, NX and NY: pixel number in x-y (matrix)\nNXY: consecutive pixel number, all the rest values are weight percents including the total.\nUnfortunatelly this script does not work with more oxides or halogens at the moment, but can be easily\nmodified to do so by the user. Just keep in mind to also modify the external functions that are being used.\n')

#Reading excel file with major and trace elements in glasses
name = input('Please type the path where the data file is located e.g. Documents/Python/Mydata.xlsx')

#data        = pd.read_excel (r'09028A_M1.xlsx',skiprows=[0])
data        = pd.read_excel (name,)
#Creating two arrays (composition and coordinates) using data from excel file using the name of the columns
df_majors   = pd.DataFrame(data, columns = ['SiO2 WT%','TiO2 WT%','Al2O3 WT%','FeO WT%','MnO WT%','MgO WT%',
                                         'CaO WT%','Na2O WT%','K2O WT%','P2O5 WT%','Total'])
df_coord    = pd.DataFrame(data, columns = ['X','Y','NX','NY'])
#transform array to a new one which can bed use to perform calculations
data_majors = df_majors.to_numpy()
data_coord  = df_coord.to_numpy()

# ..
c           = int(max(data_coord[:,2]))-1 #max num pixel on x
d           = int(max(data_coord[:,3]))-1 #max num pixel on y
#setting up empty molar fractions arrays
NK_A        = np.zeros((d,c)) #NK/A
NKC_A       = np.zeros((d,c)) #NKC/A
Mg_MgFe2    = np.zeros((d,c)) #XMg(Fe2+)
K_Na        = np.zeros((d,c)) #K/Na
#setting up empty cation fractions arrays
Mf          = np.zeros((d,c)) #M factor
#setting up empty anhydrous-based oxides wt.% arrays
SiO2        = np.zeros((d,c))
TiO2        = np.zeros((d,c))
Al2O3       = np.zeros((d,c))
FeO         = np.zeros((d,c))
MnO         = np.zeros((d,c))
MgO         = np.zeros((d,c))
CaO         = np.zeros((d,c))
Na2O        = np.zeros((d,c))
K2O         = np.zeros((d,c))
P2O5        = np.zeros((d,c))
#setting up empty Normative Differentiation Index array
DI          = np.zeros((d,c))
#setting up empty normative minerals wt.% arrays
Q           = np.zeros((d,c)) #quartz
Fsp         = np.zeros((d,c)) #albite + orthoclase
Ne          = np.zeros((d,c)) #nepheline
#setting up empty coordinate arrays (mesh-like)
Nx          = np.zeros((d,c)) #coordinates of pixels (columns)
Ny          = np.zeros((d,c)) #coordinates of pixels (rows)

# ..

data_anhf   = to_anhydrous(data_majors) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf)         #calculates anhydrous-based cation compositions mol
data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors))) #calculates oxide compositions including Fe2O3 mol
data_normf  = norm_calc(data_molf2) #calculates normative mineralogy
#reorganising calculated data based on pixel coordinates
#filling up empty arrays
for j in range (d):
    for i in range (c):
        SiO2[j,i]     = data_anhf[j*c+j+i,0]
        TiO2[j,i]     = data_anhf[j*c+j+i,1]
        Al2O3[j,i]    = data_anhf[j*c+j+i,2]
        FeO[j,i]      = data_anhf[j*c+j+i,3]
        MnO[j,i]      = data_anhf[j*c+j+i,4]
        MgO[j,i]      = data_anhf[j*c+j+i,5]
        CaO[j,i]      = data_anhf[j*c+j+i,6]
        Na2O[j,i]     = data_anhf[j*c+j+i,7]
        K2O[j,i]      = data_anhf[j*c+j+i,8]
        P2O5[j,i]     = data_anhf[j*c+j+i,9]
        Mf[j,i]       = (data_catf[j*(c+1)+i,7]+data_catf[j*(c+1)+i,8]+
                         (data_catf[j*(c+1)+i,6]*2))/(data_catf[j*(c+1)+i,0]*data_catf[j*(c+1)+i,2])
        NK_A[j,i]     = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8])/data_molf[j*(c+1)+i,2]
        NKC_A[j,i]    = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8]+
                         data_molf[j*(c+1)+i,6])/data_molf[j*(c+1)+i,2]
        Mg_MgFe2[j,i] = data_molf2[j*(c+1)+i,6]/(data_molf2[j*(c+1)+i,4]+data_molf2[j*(c+1)+i,6]+
                                                 data_molf2[j*(c+1)+i,5])
        K_Na[j,i]     = (data_anhf[j*c+j+i,8]*8302)/(data_anhf[j*c+j+i,7]*7419)
        DI[j,i]       = np.sum(data_normf[j*c+j+i,0:5])+data_normf[j*c+j+i,8]
        Q[j,i]        = data_normf[j*c+j+i,0]
        Fsp[j,i]      = data_normf[j*c+j+i,1]+data_normf[j*c+j+i,3]
        Ne[j,i]       = data_normf[j*c+j+i,4]
        Nx[j,i]       = data_coord[j*c+j+i,2]
        Ny[j,i]       = data_coord[j*c+j+i,3]

#setting up plotting scheme for SiO2, Al2O3, CaO, Na2O, K2O, K/Na
cm       = 1/2.54
plt.figure(figsize = (6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(Fsp,cmap='Spectral',vmin=0, vmax=100)
axs[0,0].set_title("SiO$_{2}$ wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Q,cmap='Spectral',vmin=0, vmax=40)
axs[0,1].set_title("Al$_{2}$O$_{3}$ wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(CaO,cmap='Spectral',vmin=0, vmax=60)
axs[0,2].set_title("CaO wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(Na2O,cmap='Spectral',vmin=0, vmax=15)
axs[1,0].set_title("Na$_{2}$O wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(K2O,cmap='Spectral',vmin=0, vmax=15)
axs[1,1].set_title("K$_{2}$O wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(K_Na,cmap='Spectral',vmin=0, vmax=5)
axs[1,2].set_title("K/Na")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig('09028A_M1a.pdf',dpi=600, transparent=True, bbox_inches='tight')
