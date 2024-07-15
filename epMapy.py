#Importing libraries a$nd changing plotting parameters

import numpy as np #math functions as pi,e,random,arrays, etc 
import matplotlib.pyplot as plt #default plotting library
import pandas as pd #for printing nice tables (pd.DataFrame)
import seaborn as sbn #extra version of plt
import time 
plt.rcParams['pdf.fonttype'] = 42 #to use Type 42 (a.k.a. TrueType) fonts for PDF files
plt.rcParams['ps.fonttype']  = 42 #to use Type 42 (a.k.a. TrueType) fonts for PostScript files
plt.rcParams['svg.fonttype'] = 'none' #for neither embedding the font nor rendering the text as path
plt.rcParams.update({'font.size': 7}) #defined font size
plt.rcParams['figure.dpi']   = 200 #increases the quality of images
import mpltern #to plot ternary diagrams
from matplotlib_scalebar.scalebar import ScaleBar #to add a scale bar to the maps
from functionsmapping import to_anhydrous, to_mol, to_cat, norm_calc, add_fe2o3, clean_data #importing external functions

print('\n\nThis script allows plotting EPMA quantitative maps and do calculations with the data within.\nPlease be aware that this script is designed for the output file obtained with CalcImage (J. Donovan).\nThe input file for this script is a .xslx file, so you should convert .DAT to .xslx before\nat least the following column names (case and character sensitive):\nX,Y,NX,NY,NXY,SiO2 WT%,TiO2 WT%,Al2O3 WT%,FeO WT%,MnO WT%,MgO WT%,CaO WT%,Na2O WT%,K2O WT%,P2O5 WT%,Total\nX and Y: x-y coordinates EPMA position, NX and NY: pixel number in x-y (matrix)\nNXY: consecutive pixel number, all the rest values are weight percents including the total.\nUnfortunatelly this script does not work with more oxides or halogens at the moment, but can be easily\nmodified to do so by the user. Just keep in mind to also modify the external functions that are being used.\n')

#Reading excel file with major and trace elements in glasses
name = input('Please type the path where the data file is located:   ')
filename = input('Please type the file name for the plots:   ')
st = time.time()
#data        = pd.read_excel (r'09028A_M1.xlsx',skiprows=[0])
data        = pd.read_excel (name,)
#Creating two arrays (composition and coordinates) using data from excel file using the name of the columns
df_majors   = pd.DataFrame(data, columns = ['SiO2 WT%','TiO2 WT%','Al2O3 WT%','FeO WT%','MnO WT%','MgO WT%',
                                         'CaO WT%','Na2O WT%','K2O WT%','P2O5 WT%','Total'])
df_coord    = pd.DataFrame(data, columns = ['X','Y','NX','NY'])
#transform array to a new one which can bed use to perform calculations
data_majors = df_majors.to_numpy()
data_coord  = df_coord.to_numpy()

#setting up empty arrays
c           = int(max(data_coord[:,2]))-1 #max num pixel on x
d           = int(max(data_coord[:,3]))-1 #max num pixel on y
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
data_majors = clean_data(data_majors)   #cleaning holes and cracks that retrive very low totals
data_anhf   = to_anhydrous(data_majors) #calculates anhydrous-base oxide compositions wt.%
data_molf   = to_mol(data_anhf)         #calculates anhydrous-based oxide compositions mol
data_catf   = to_cat(data_molf)         #calculates anhydrous-based cation compositions mol
data_molf2  = to_mol(to_anhydrous(add_fe2o3(data_majors))) #calculates oxide compositions including Fe2O3 mol
data_normf  = norm_calc(data_molf2) #calculates normative mineralogy

#filtered_df = df[df['emails'].str.contains('GMAIL', case=False)]
#reorganising calculated data based on pixel coordinates
#filling up empty arrays
for j in range (d):
    for i in range (c):
        if data_majors[j*c+j+i,-1]>0:
            SiO2p[j,i],TiO2p[j,i],Al2O3p[j,i],FeOp[j,i],MnOp[j,i],MgOp[j,i],CaOp[j,i],Na2Op[j,i],K2Op[j,i],P2O5p[j,i],totp[j,i]     = [data_majors[j*c+j+i,q] for q in range (11)]
            SiO2[j,i],TiO2[j,i],Al2O3[j,i],FeO[j,i],MnO[j,i],MgO[j,i],CaO[j,i],Na2O[j,i],K2O[j,i],P2O5[j,i],tot[j,i]                = [data_anhf[j*c+j+i,q] for q in range (11)]
            Q[j,i],ort[j,i],lc[j,i],ab[j,i],ne[j,i],an[j,i],C[j,i],ac[j,i],ns[j,i],di[j,i],ol[j,i],hy[j,i],mt[j,i],ilm[j,i],ap[j,i] = [data_normf[j*c+j+i,q] for q in range (15)]
            Mf[j,i]        = (data_catf[j*(c+1)+i,7]+data_catf[j*(c+1)+i,8]+
                             (data_catf[j*(c+1)+i,6]*2))/(data_catf[j*(c+1)+i,0]*data_catf[j*(c+1)+i,2])
            NK_A[j,i]      = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8])/data_molf[j*(c+1)+i,2]
            NKC_A[j,i]     = (data_molf[j*(c+1)+i,7]+data_molf[j*(c+1)+i,8]+
                             data_molf[j*(c+1)+i,6])/data_molf[j*(c+1)+i,2]
            Mg_MgFe2[j,i]  = data_molf2[j*(c+1)+i,6]/(data_molf2[j*(c+1)+i,4]+data_molf2[j*(c+1)+i,6]+
                             data_molf2[j*(c+1)+i,5])
            K_Na[j,i]      = (data_anhf[j*c+j+i,8]*8302)/(data_anhf[j*c+j+i,7]*7419)
            DI[j,i]        = np.sum(data_normf[j*c+j+i,0:5])+data_normf[j*c+j+i,8]
            Fsp[j,i]       = ort[j,i]+ab[j,i]+an[j,i]
            px[j,i]        = ac[j,i]+di[j,i]+hy[j,i]
            ox[j,i]        = mt[j,i]+ilm[j,i]
            Nx[j,i]        = data_coord[j*c+j+i,2]
            Ny[j,i]        = data_coord[j*c+j+i,3]
        else:
            SiO2p[j,i],TiO2p[j,i],Al2O3p[j,i],FeOp[j,i],MnOp[j,i],MgOp[j,i],CaOp[j,i],Na2Op[j,i],K2Op[j,i],P2O5p[j,i]               = [300 for q in range (10)]
            SiO2[j,i],TiO2[j,i],Al2O3[j,i],FeO[j,i],MnO[j,i],MgO[j,i],CaO[j,i],Na2O[j,i],K2O[j,i],P2O5[j,i]                         = [300 for q in range (10)]
            totp[j,i],tot[j,i],Mf[j,i],K_Na[j,i]                                                                                    = [-1 for q in range (4)]
            NK_A[j,i],NKC_A[j,i],Mg_MgFe2[j,i],DI[j,i]                                                                              = [300 for q in range (4)]
            Q[j,i],ort[j,i],lc[j,i],ab[j,i],ne[j,i],an[j,i],C[j,i],ac[j,i],ns[j,i],di[j,i],ol[j,i],hy[j,i],mt[j,i],ilm[j,i],ap[j,i] = [300 for q in range (15)]
            Nx[j,i]        = data_coord[j*c+j+i,2]
            Ny[j,i]        = data_coord[j*c+j+i,3]


for j in range (d):
    for i in range (c):
        if (Q[j,i] <20) & (Fsp[j,i]>40) & (px[j,i]>10):
            mrg[j,i]=5
        elif (Q[j,i] >80):
            mrg[j,i]=1
        elif (ox[j,i] >60):
            mrg[j,i]=3
        elif (px[j,i] >70):
            mrg[j,i]=4
        elif (ap[j,i] >80):
            mrg[j,i]=5

#setting up plotting scheme for SiO2, Al2O3, CaO, Na2O, K2O, K/Na
cm       = 1/2.54
plt.figure(figsize = (6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(Q,cmap='Spectral_r',vmin=0, vmax=100)
a0.cmap.set_over('black')
axs[0,0].set_title("Q wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Fsp,cmap='Spectral_r',vmin=0, vmax=100)
a1.cmap.set_over('black')
axs[0,1].set_title("ab + or + an wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(px,cmap='Spectral_r',vmin=0, vmax=100)
a2.cmap.set_over('black')
axs[0,2].set_title("px wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(ol,cmap='Spectral_r',vmin=0, vmax=100)
a3.cmap.set_over('black')
axs[1,0].set_title("ol wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(ox,cmap='Spectral_r',vmin=0, vmax=100)
a4.cmap.set_over('black')
axs[1,1].set_title("mt + ilm wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(ap,cmap='Spectral_r',vmin=0, vmax=100)
a5.cmap.set_over('black')
axs[1,2].set_title("ap  wt.%")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig(filename+'_normative.pdf',dpi=600, transparent=True, bbox_inches='tight')

#setting up plotting scheme for TiO2, P2O5, FeO, MgO, differentiation index, M-factor
plt.figure(figsize=(6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(SiO2,cmap='Spectral_r',vmin=0, vmax=100)
a0.cmap.set_over('black')
axs[0,0].set_title("SiO$_{2}$ wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Al2O3,cmap='Spectral_r',vmin=0, vmax=45)
a1.cmap.set_over('black')
axs[0,1].set_title("Al$_{2}$O$_{3}$ wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(FeO,cmap='Spectral_r',vmin=0, vmax=100)
a2.cmap.set_over('black')
axs[0,2].set_title("FeO wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(MgO,cmap='Spectral_r',vmin=0, vmax=45)
a3.cmap.set_over('black')
axs[1,0].set_title("MgO wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(Na2O,cmap='Spectral_r',vmin=0, vmax=15)
a4.cmap.set_over('black')
axs[1,1].set_title("Na$_{2}$O wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(P2O5,cmap='Spectral_r',vmin=0, vmax=60)
a5.cmap.set_over('black')
axs[1,2].set_title("P$_{2}$O$_{5}$")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig(filename+'_oxides_anh.pdf',dpi=600, transparent=True, bbox_inches='tight')

plt.figure(figsize=(6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(SiO2p,cmap='Spectral_r',vmin=0, vmax=100)
a0.cmap.set_over('black')
axs[0,0].set_title("SiO$_{2}$ wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Al2O3p,cmap='Spectral_r',vmin=0, vmax=45)
a1.cmap.set_over('black')
axs[0,1].set_title("Al$_{2}$O$_{3}$ wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(FeOp,cmap='Spectral_r',vmin=0, vmax=100)
a2.cmap.set_over('black')
axs[0,2].set_title("FeO wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(MgOp,cmap='Spectral_r',vmin=0, vmax=45)
a3.cmap.set_over('black')
axs[1,0].set_title("MgO wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(Na2Op,cmap='Spectral_r',vmin=0, vmax=15)
a4.cmap.set_over('black')
axs[1,1].set_title("Na$_{2}$O wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(totp,cmap='Spectral_r',vmin=1, vmax=100)
a5.cmap.set_under('black')
axs[1,2].set_title("Total wt.%")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig(filename+'_oxides.pdf',dpi=600, transparent=True, bbox_inches='tight')

#setting up plotting scheme for NK/A, NKC/A, XMg, normative Q, normative ab+or, normative Ne
cm       = 1/2.54
plt.figure(figsize=(6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(NK_A,cmap='Spectral_r',vmin=0, vmax=2)
a0.cmap.set_over('black')
axs[0,0].set_title("NK/A")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(mrg,cmap='tab10')
a1.cmap.set_over('black')
axs[0,1].set_title("PhaseMap")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(Mf,cmap='Spectral_r',vmin=0, vmax=5)
a2.cmap.set_under('black')
axs[0,2].set_title("M$_{f}$")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(Mg_MgFe2,cmap='Spectral_r',vmin=0, vmax=1)
a3.cmap.set_over('black')
axs[1,0].set_title("XMg$_{Fe2+}$")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(K_Na,cmap='Spectral_r',vmin=0, vmax=5)
a4.cmap.set_under('black')
axs[1,1].set_title("K/Na")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(DI,cmap='Spectral_r',vmin=0, vmax=100)
a5.cmap.set_over('black')
axs[1,2].set_title("Di")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig(filename+'_indexes.pdf',dpi=600, transparent=True, bbox_inches='tight')


#Reshaping matrix into 1D-arrays and deleting nan values
K_Na_new = K_Na.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
Mf_new = Mf.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
SiO2_new = SiO2.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
FeO_new = FeO.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
CaO_new = CaO.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
P2O5_new = P2O5.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
Nx_new = Nx.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
Ny_new = Ny.reshape([1, c*d])[~np.isnan(K_Na.reshape([1, c*d]))]
#Filtering K_Na data to only consider rhyolitic glass
#sbn.kdeplot(K_Na_new,bw_adjust=.4)
K_Na_new2=K_Na_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
Mf_new2=Mf_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
SiO2_new2=SiO2_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
FeO_new2=FeO_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
CaO_new2=CaO_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
P2O5_new2=P2O5_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
Nx_new2=Nx_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
Ny_new2=Ny_new[(SiO2_new > 72) & (SiO2_new < 84) & (FeO_new < 2.5) & (CaO_new < 4) &
                   (K_Na_new < 6) & (K_Na_new > 0) & (P2O5_new < 1) & (P2O5_new > 0)]
#setting up plotting scheme for NK/A, NKC/A, XMg, normative Q, normative ab+or, normative Ne
cm       = 1/2.54
plt.figure(figsize=(6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
sbn.histplot(K_Na_new2,stat='percent', kde=True, bins=30,ax=axs[0,0])
axs[0,0].set_xlabel("K/N")
axs[0,0].set_ylabel(None)
sbn.histplot(SiO2_new2,stat='percent', kde=True,bins=30, ax=axs[0,1])
axs[0,1].set_xlabel("SiO$_{2}$ wt.%")
axs[0,1].set_ylabel(None)
sbn.histplot(FeO_new2,stat='percent', kde=True, bins=30,ax=axs[0,2])
axs[0,2].set_xlabel("FeO wt.%")
axs[0,2].set_ylabel(None)
sbn.histplot(CaO_new2,stat='percent', kde=True,bins=30, ax=axs[1,0])
axs[1,0].set_xlabel("CaO wt.%")
axs[1,0].set_ylabel(None)
sbn.histplot(P2O5_new2,stat='percent', kde=True, kde_kws={'bw_adjust':3},bins=30, ax=axs[1,1])
axs[1,1].set_xlabel("P$_{2}$O$_{5}$ wt.%")
axs[1,1].set_ylabel(None)
sbn.histplot(Mf_new2,stat='percent', kde=True,bins=30, ax=axs[1,2])
axs[1,2].set_xlabel("M$_{f}$")
axs[1,2].set_ylabel(None)
plt.tight_layout()
plt.savefig(filename+'_rastered1.pdf',dpi=600, transparent=True, bbox_inches='tight')

Nx_new3=Nx_new2[(Mf_new2 < 1.2)&(K_Na_new2 < 2.5)&(K_Na_new2 > 0.5)]
Ny_new3=Ny_new2[(Mf_new2 < 1.2)&(K_Na_new2 < 2.5)&(K_Na_new2 > 0.5)]
#setting up plotting scheme for SiO2, Al2O3, CaO, Na2O, K2O, K/Na
cm       = 1/2.54
plt.figure(figsize = (6.5*cm,6.5*cm))
fig, axs = plt.subplots(2, 3)
a0       = axs[0,0].imshow(SiO2,cmap='Spectral_r',vmin=0, vmax=100)
a0.cmap.set_over('black')
axs[0,0].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[0,0].set_title("SiO$_{2}$ wt.%")
axs[0,0].axis("off")
fig.colorbar(a0)
a1       = axs[0,1].imshow(Al2O3,cmap='Spectral_r',vmin=0, vmax=40)
a1.cmap.set_over('black')
axs[0,1].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[0,1].set_title("Al$_{2}$O$_{3}$ wt.%")
axs[0,1].axis("off")
fig.colorbar(a1)
a2       = axs[0,2].imshow(CaO,cmap='Spectral_r',vmin=0, vmax=60)
a2.cmap.set_over('black')
axs[0,2].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[0,2].set_title("CaO wt.%")
axs[0,2].axis("off")
fig.colorbar(a2)
a3       = axs[1,0].imshow(Na2O,cmap='Spectral_r',vmin=0, vmax=15)
a3.cmap.set_over('black')
axs[1,0].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[1,0].set_title("Na$_{2}$O wt.%")
axs[1,0].axis("off")
fig.colorbar(a3)
a4       = axs[1,1].imshow(K2O,cmap='Spectral_r',vmin=0, vmax=15)
a4.cmap.set_over('black')
axs[1,1].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[1,1].set_title("K$_{2}$O wt.%")
axs[1,1].axis("off")
fig.colorbar(a4)
a5       = axs[1,2].imshow(K_Na,cmap='Spectral_r',vmin=0, vmax=5)
a5.cmap.set_under('black')
axs[1,2].scatter(Nx_new3,Ny_new3,c='k',marker='o',edgecolors="white",linewidths=0.25,s=2)
axs[1,2].set_title("K/Na")
axs[1,2].axis("off")
fig.colorbar(a5)
scalebar = ScaleBar(0.000005) # 1 pixel = 5 µm
plt.gca().add_artist(scalebar)
plt.tight_layout()
plt.savefig(filename+'_rastered2.pdf',dpi=600, transparent=True, bbox_inches='tight')
et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')