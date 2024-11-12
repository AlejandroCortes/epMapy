#importing libraries for standalone functions
import numpy as np
import pandas as pd
import re

def clean_data(X, oxide_column_indices):
    a = len(X[:, 0])  # number of data points
    b = len(X[0, :])  # number of oxides + total
    clean = np.copy(X)  # copies the input to avoid overwriting

    # Accessing the FeO, CaO, and Total columns using the indices
    feo_index = oxide_column_indices.get('FeO', None)
    cao_index = oxide_column_indices.get('CaO', None)
    total_index = oxide_column_indices.get('Total', None)

    for i in range(a):
        # Retrieve the values dynamically using column indices
        if (clean[i, total_index] < 80) & ((clean[i, feo_index] < 40) | (clean[i, cao_index] < 50)):
            clean[i, :] = clean[i, :] * 0  # flagging as a bad analysis or crackj
    
    return clean

#This function takes major oxide compositions and normalise them to 100%
def to_anhydrous(X, oxide_column_indices):
    a   = len(X[:,0]) #num of data points
    b   = len(X[0,:]) #num of oxides + total
    anh = np.copy(X)  #copies the input to avoid overwriting
    total_index = oxide_column_indices.get('Total', None)

    for i in range (a):
        for j in range (b):
            if (anh[i, total_index] > 0):
                anh[i,j] = anh[i,j]*100/anh[i,total_index] #normalising data to 100%
    return anh

#This function takes major oxide compositions in wt.% and converts them to mol
def to_mol(X, oxide_column_indices):
    a = len(X[:, 0])  # number of pixels (rows)
    b = len(X[0, :])  # number of oxides + total (columns)
    total_index = oxide_column_indices.get('Total', None)

    # Molecular weights with FeO and Fe2O3 calculated
    mw = {
        "SiO2": 60.080, "TiO2": 79.870, "Cr2O3": 152.000, "Al2O3": 101.960, "FeO": 71.840, 
        "MnO": 70.940, "NiO": 74.690, "MgO": 40.300, "CaO": 56.080, "Na2O": 61.980, 
        "K2O": 94.200, "P2O5": 283.880, "F": 18.998, "Cl": 35.450, "SO3": 80.060, 
        "BaO": 153.330, "SrO": 103.620, "Fe2O3": 141.94
    }

    mol = np.copy(X)  # Copy the input to avoid overwriting

    # Iterate through the data to convert oxides to mol
    for i in range(a):  # Iterate over each row (data point)
        for col_name, col_index in oxide_column_indices.items():  # Iterate over oxide columns
            if mol[i, total_index] > 0:  # Skip if the value is zero or negative
                if col_name == "Total":
                    mol[i, col_index] = 100  # Set "Total" to 100
                else:
                    mol[i, col_index] = mol[i, col_index] / mw.get(col_name, 1)  # Convert to mol using molecular weight
    return mol

#This function takes major oxide compositions in mol and converts them to cation fractions
def to_cat(X, oxide_column_indices):
    a = len(X[:, 0])  # number of pixels
    b = len(X[0, :])  # number of oxides + total
    total_index = oxide_column_indices.get('Total', None)
    
    ct = {
        "SiO2": 1, "TiO2": 1, "Cr2O3": 2, "Al2O3": 2, "FeO": 1,
        "MnO": 1, "NiO": 1, "MgO": 1, "CaO": 1, "Na2O": 2,
        "K2O": 2, "P2O5": 2, "F": 1, "Cl": 1, "SO3": 1,
        "BaO": 1, "SrO": 1, "Fe2O3": 2
    }
    
    cat = np.copy(X)  # copies the input to avoid overwriting
    
    for i in range(a):
        # Calculate cation fractions for each oxide column
        if cat[i, total_index] > 0:  # Skip if the total value is zero or negative
            for col_name, col_index in oxide_column_indices.items():  # Iterate over oxide columns
                if col_name != "Total":
                    cat[i, col_index] = cat[i, col_index] * ct.get(col_name, 1)  # Calculate cation fractions
        
            # Sum all columns (excluding 'Total' column) and store in 'Total' column
            oxide_columns = [index for col_name, index in oxide_column_indices.items() if col_name != "Total"]  # Exclude "Total"
            sum_of_oxides = np.sum(cat[i, oxide_columns])  # sum all columns except 'Total'
            cat[i, total_index] = sum_of_oxides  # Store the sum in the 'Total' column
            for col_name, col_index in oxide_column_indices.items():
                if col_name != "Total":
                    cat[i, col_index] = cat[i, col_index] / cat[i, total_index] 
                
    return cat

def add_fe2o3(X, oxide_column_indices):
    Fe2_Fetot = 0.7  # Ratio Fe2+/Fe3+ for NNO
    X2 = np.copy(X)  # Copy the input to avoid overwriting

    # Check if 'Fe2O3' column already exists in oxide_column_indices
    if "Fe2O3" not in oxide_column_indices:
        # Calculate Fe2O3 from FeOtot and Fe2+/Fetot ratio
        feo_index = oxide_column_indices.get('FeO', None)
        Fe2O3 = X[:, feo_index] * (1 - Fe2_Fetot) * 1.11134
        
        # Add the Fe2O3 data to the array
        X2 = np.append(X2, np.reshape(Fe2O3, (len(Fe2O3), 1)), axis=1)
        
        # Update FeO with the Fe2+/Fetot ratio
        X2[:, feo_index] = X2[:, feo_index] * Fe2_Fetot
    
        # Update the oxide_column_indices to include Fe2O3
        oxide_column_indices["Fe2O3"] = len(oxide_column_indices)  # Add Fe2O3 at the next available index
    else:
        # If Fe2O3 already exists, just return the original array without changes
        X2 = np.copy(X)  # No changes to the input array
        print("Fe2O3 column already exists, no update performed.")
    
    return X2, oxide_column_indices

def norm_calc(X, oxide_column_indices):
    Z = np.copy(X)  # Copy the input to avoid overwriting
    a = len(Z[:, 0])  # Number of pixels
    oots = np.zeros((a, 23))  # Array to store calculated values
    fps=np.zeros((a,8)) #First pass silica
    sps=np.zeros((a,11)) #Second pass silica
    fs=np.zeros((a,7)) #final silica
    n_wt=np.zeros((a,16)) #recalculate to normative mineral in wt%

    # Define the oxides we're interested in and map them to indices in oxide_column_indices
    oxide_names = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O','P2O5', 'Total']
    oots_names = [
        "ap_Ca", "ap_P", "ilm_Fe_II", "Ilm_Ti", "or_lc_Al", "or_lc_K", "ab_ne_Al", "ab_ne_Na", 
        "ac_Na", "ac_FeIII", "ns_Na", "mt_Fe_II", "mt_Fe_III", "an_Al_prelim", "an_Al_final", 
        "C_Al", "an_Ca", "di_Ca", "di_Mg", "di_FeII", "ol_hy_Fe_II", "ol_hy_Mg", "Mg_num"
        ]
    fps_names = ["lc", "ne", "ac", "ns", "an", "di", "ol", "residual_si"]
    sps_names = [
        "or_K", "lc_K", "residual_si2", "ab_Na", "residual_si3", "prelim_ol_Mg_Fe", "prelim_hy_Mg_Fe", 
        "actual_ol_Mg_Fe", "actual_hy_Mg_Fe", "residual_siQ"
        ]
    fs_names = ["or", "lc", "ab", "ne", "ol", "hy", "Q"]
    nwt_names = ["Q", "or", "lc", "ab", "ne", "an", "C", "ac", "ns", "di", "ol", "hy", "mt", "ilm", "ap","total"]

    mw = {
        "SiO2": 60.080, "TiO2": 79.870, "Cr2O3": 152.000, "Al2O3": 101.960, "FeO": 71.840, 
        "MnO": 70.940, "NiO": 74.690, "MgO": 40.300, "CaO": 56.080, "Na2O": 61.980, 
        "K2O": 94.200, "P2O5": 283.880, "F": 18.998, "Cl": 35.450, "SO3": 80.060, 
        "BaO": 153.330, "SrO": 103.620, "Fe2O3": 141.94
    }

    oxide_indices = {oxide: oxide_column_indices.get(oxide, None) for oxide in oxide_names}
    oots_index = {name: index for index, name in enumerate(oots_names)}
    fps_index = {name: index for index, name in enumerate(fps_names)}
    sps_index = {name: index for index, name in enumerate(sps_names)}
    fs_index = {name: index for index, name in enumerate(fs_names)}
    nwt_index = {name: index for index, name in enumerate(nwt_names)}

    condA = Z[:, oxide_indices['Total']] > 0

    oots[:, oots_index["ap_Ca"]] = Z[:, oxide_indices['P2O5']] * 10 / 3  # ap Ca
    oots[:, oots_index["ap_P"]] = Z[:, oxide_indices['P2O5']]  # ap P
    oots[:, oots_index["ilm_Fe_II"]] = Z[:, oxide_indices['TiO2']]  # ilm Fe(II)
    oots[:, oots_index["Ilm_Ti"]] = Z[:, oxide_indices['TiO2']]  # ilm Ti
    oots[:, oots_index["or_lc_Al"]] = Z[:, oxide_indices['K2O']]  # or/lc Al
    oots[:, oots_index["or_lc_K"]] = Z[:, oxide_indices['K2O']]  # or/lc K

    oots[:, oots_index["ab_ne_Al"]] = np.where(condA & (Z[:, oxide_indices['Na2O']] > (Z[:, oxide_indices['Al2O3']] - oots[:, oots_index["or_lc_Al"]])),
                                               Z[:, oxide_indices['Al2O3']] - oots[:, oots_index["or_lc_Al"]],
                                               Z[:, oxide_indices['Na2O']]) 
    oots[:, oots_index["ab_ne_Na"]] = np.where(condA,
                                                oots[:, oots_index["ab_ne_Al"]])
    oots[:, oots_index["ac_FeIII"]] = np.where(condA & ((Z[:, oxide_indices['Na2O']] - oots[:, oots_index["ab_ne_Na"]]) < Z[:, oxide_indices['Fe2O3']]),
                                               Z[:, oxide_indices['Na2O']] - oots[:, oots_index["ab_ne_Na"]],
                                               Z[:, oxide_indices['Fe2O3']])
    oots[:, oots_index["ac_Na"]] = np.where(condA, oots[:, oots_index["ac_FeIII"]])  # ac Na
    oots[:, oots_index["ns_Na"]] = np.where(condA, Z[:, oxide_indices['Na2O']] - oots[:, oots_index["ab_ne_Na"]] - oots[:, oots_index["ac_Na"]])  # ns Na
    oots[:, oots_index["mt_Fe_III"]] = np.where(condA, Z[:, oxide_indices['Fe2O3']] - oots[:, oots_index["ac_FeIII"]])  # mt Fe(III)
    oots[:, oots_index["mt_Fe_II"]] = np.where(condA, oots[:, oots_index["mt_Fe_III"]])
    oots[:, oots_index["an_Al_prelim"]] = np.where(condA & (oots[:, oots_index["ac_Na"]] > 0.000001),
                                                   0,
                                                   Z[:, oxide_indices['Al2O3']] - oots[:, oots_index["or_lc_Al"]] - oots[:, oots_index["ab_ne_Al"]])
    oots[:, oots_index["an_Al_final"]] = np.where(condA & (oots[:, oots_index["an_Al_prelim"]] > (Z[:, oxide_indices['CaO']] - oots[:, oots_index["ap_Ca"]])),
                                                  Z[:, oxide_indices['CaO']] - oots[:, oots_index["ap_Ca"]],
                                                  oots[:, oots_index["an_Al_prelim"]])                                               
    oots[:, oots_index["C_Al"]] = np.where(condA &(oots[:, oots_index["an_Al_prelim"]] == oots[:, oots_index["an_Al_final"]]),
                                           0,
                                           oots[:, oots_index["an_Al_prelim"]] - oots[:, oots_index["an_Al_final"]])
    oots[:, oots_index["an_Ca"]] = np.where(condA, oots[:, oots_index["an_Al_final"]])  # an Ca
    oots[:, oots_index["di_Ca"]] = np.where(condA, Z[:, oxide_indices['CaO']] - oots[:, oots_index["ap_Ca"]] - oots[:, oots_index["an_Ca"]])  # di Ca
    oots[:, oots_index["Mg_num"]] = np.where(condA, Z[:, oxide_indices['MgO']] / (Z[:, oxide_indices['MgO']] + Z[:, oxide_indices['FeO']] + 
                                                                                                             Z[:, oxide_indices['MgO']] - oots[:, oots_index["ilm_Fe_II"]] - 
                                                                                                             oots[:, oots_index["mt_Fe_II"]]))  # XMg
    oots[:, oots_index["di_Mg"]] = np.where(condA, oots[:, oots_index["di_Ca"]] * oots[:, oots_index["Mg_num"]])  # di Mg
    oots[:, oots_index["di_FeII"]] = np.where(condA, oots[:, oots_index["di_Ca"]] * (1 - oots[:, oots_index["Mg_num"]]))  # di FeII
    oots[:, oots_index["ol_hy_Fe_II"]] = np.where(condA, Z[:, oxide_indices['FeO']] + Z[:, oxide_indices['MnO']] - oots[:, oots_index["ilm_Fe_II"]] - 
                                                  oots[:, oots_index["mt_Fe_II"]] - oots[:, oots_index["di_FeII"]])  # ol/hy Fe(II)
    oots[:, oots_index["ol_hy_Mg"]] = np.where(condA,Z[:, oxide_indices['MgO']] - oots[:, oots_index["di_Mg"]])  # ol/hy Mg
    
    fps[:, fps_index["lc"]]   = np.where(condA,4*oots[:, oots_index["or_lc_K"]])                               #lc
    fps[:, fps_index["ne"]]   = np.where(condA,2*oots[:, oots_index["ab_ne_Na"]])                               #ne
    fps[:, fps_index["ac"]]   = np.where(condA,4*oots[:, oots_index["ac_Na"]])                                        #ac
    fps[:, fps_index["ns"]]   = np.where(condA,oots[:, oots_index["ns_Na"]])                                         #ns
    fps[:, fps_index["an"]]   = np.where(condA,2*oots[:, oots_index["an_Ca"]])                                       #an
    fps[:, fps_index["di"]]   = np.where(condA,2*oots[:, oots_index["di_Ca"]])                                       #di
    fps[:, fps_index["ol"]]   = np.where(condA,(oots[:, oots_index["ol_hy_Fe_II"]]+ oots[:, oots_index["ol_hy_Mg"]])/2)  
    
    fps_columns = [index for col_name, index in fps_index.items() if col_name != "residual_si"]
    sum_of_fps = np.where(condA,np.sum(fps[:, fps_columns], axis=1))  # sum all columns except "residual_si"
    fps[:, fps_index["residual_si"]]   = np.where(condA,Z[:, oxide_indices["SiO2"]]-sum_of_fps) #residual Si
    
    sps[:, sps_index["or_K"]] = np.where(condA & ((oots[:, oots_index["or_lc_K"]]*2) > fps[:, fps_index["residual_si"]]),
                                        (fps[:, fps_index["lc"]]+(fps[:,fps_index["residual_si"]]-(oots[:, oots_index["or_lc_K"]]*2)))/2, #or K
                                         oots[:, oots_index["or_lc_K"]])                             #or K 
    sps[:, sps_index["lc_K"]]  = np.where(condA,oots[:, oots_index["or_lc_K"]]-sps[:, sps_index["or_K"]])                         #lc K
    sps[:, sps_index["residual_si2"]]  = np.where(condA,fps[:, fps_index["residual_si"]]-(2*sps[:, sps_index["or_K"]]))                      #residual Si
    sps[:, sps_index["ab_Na"]] = np.where(condA & ((oots[:, oots_index["ab_ne_Na"]]*4) > sps[:, sps_index["residual_si2"]]),
                                          (sps[:, sps_index["residual_si2"]]+fps[:, fps_index["ne"]]-(2*oots[:, oots_index["ab_ne_Na"]]))/4,
                                          oots[:, oots_index["ab_ne_Na"]])
    sps[:, sps_index["ne_Na"]]  = np.where(condA,oots[:, oots_index["ab_ne_Na"]]-sps[:, sps_index["ab_Na"]])                         #ne Na
    sps[:, sps_index["residual_si3"]]  = np.where(condA,sps[:, sps_index["residual_si2"]]-(4*sps[:, sps_index["ab_Na"]]))                      #residual Si
    sps[:, sps_index["prelim_ol_Mg_Fe"]] = np.where(condA & (sps[:, sps_index["ne_Na"]] > 0.0000001),
                                                    oots[:, oots_index["ol_hy_Fe_II"]]+oots[:, oots_index["ol_hy_Mg"]],
                                                    oots[:, oots_index["ol_hy_Fe_II"]]+oots[:, oots_index["ol_hy_Mg"]]-(2*sps[:, sps_index["prelim_ol_Mg_Fe"]]))
    sps[:, sps_index["prelim_hy_Mg_Fe"]]  = np.where(condA,oots[:, oots_index["ol_hy_Fe_II"]]+oots[:, oots_index["ol_hy_Mg"]]-sps[:, sps_index["prelim_ol_Mg_Fe"]])             #prelim hy Mg+Fe
    condB = sps[:, sps_index["prelim_ol_Mg_Fe"]] < 0
    sps[:, sps_index["actual_ol_Mg_Fe"]] = np.where(condA & condB,
                                                    0,
                                                    sps[:, sps_index["prelim_ol_Mg_Fe"]])
    sps[:, sps_index["actual_hy_Mg_Fe"]] = np.where(condA & condB,
                                                    oots[:, oots_index["ol_hy_Fe_II"]]+oots[:, oots_index["ol_hy_Mg"]],
                                                    sps[:, sps_index["prelim_hy_Mg_Fe"]])
    sps[:, sps_index["residual_siQ"]] = np.where(condA & condB,
                                                 sps[:, sps_index["residual_si3"]]-((oots[:, oots_index["ol_hy_Fe_II"]]+oots[:, oots_index["ol_hy_Mg"]])/2),
                                                 0)
                                                 
    fs[:, fs_index["or"]]  = np.where(condA,6*sps[:, sps_index["or_K"]])   #or
    fs[:, fs_index["lc"]]  = np.where(condA,4*sps[:, sps_index["lc_K"]])   #lc
    fs[:, fs_index["ab"]]  = np.where(condA,6*sps[:, sps_index["ab_Na"]] )  #ab
    fs[:, fs_index["ne"]]  = np.where(condA,2*sps[:, sps_index["ne_Na"]] )  #ne
    fs[:, fs_index["ol"]]  = np.where(condA,0.5*sps[:, sps_index["actual_ol_Mg_Fe"]]) #ol
    fs[:, fs_index["hy"]]  = np.where(condA,sps[:, sps_index["actual_hy_Mg_Fe"]] )    #hy
    fs[:, fs_index["Q"]]  = np.where(condA,sps[:,sps_index["residual_siQ"]])    #Q
           
    n_wt[:, nwt_index["Q"]]  = np.where(condA,fs[:, fs_index["Q"]]*60.08)                                                 #Q
    n_wt[:, nwt_index["or"]]  = np.where(condA,sps[:, sps_index["or_K"]]*556.64)                                               #or
    n_wt[:, nwt_index["lc"]]  = np.where(condA,sps[:, sps_index["lc_K"]]*436.48)                                               #lc
    n_wt[:, nwt_index["ab"]]  = np.where(condA,sps[:, sps_index["ab_Na"]]*524.42)                                               #ab
    n_wt[:, nwt_index["ne"]]  = np.where(condA,sps[:, sps_index["ne_Na"]]*284.1)                                                #ne
    n_wt[:, nwt_index["an"]]  = np.where(condA,oots[:, oots_index["an_Ca"]]*278.2 )                                             #an
    n_wt[:, nwt_index["C"]]  = np.where(condA,oots[:, oots_index["C_Al"]]*101.96 )                                            #C
    n_wt[:, nwt_index["ac"]]  = np.where(condA,oots[:, oots_index["ac_Na"]]*462  )                                               #ac
    n_wt[:, nwt_index["ns"]]  = np.where(condA,oots[:, oots_index["ns_Na"]]*122.06  )                                           #ns
    n_wt[:, nwt_index["di"]]  = np.where(condA,(oots[:, oots_index["di_Ca"]]*116.16)+(oots[:, oots_index["di_Mg"]]*100.38)+(oots[:, oots_index["di_FeII"]]*131.39))   #di
    n_wt[:, nwt_index["ol"]] = np.where(condA,(sps[:, sps_index["actual_ol_Mg_Fe"]]*70.34*oots[:, oots_index["Mg_num"]])+(sps[:, sps_index["actual_ol_Mg_Fe"]]*(1-oots[:, oots_index["Mg_num"]])*101.89))  #ol
    n_wt[:, nwt_index["hy"]] = np.where(condA,(sps[:, sps_index["actual_hy_Mg_Fe"]]*100.38*oots[:, oots_index["Mg_num"]])+(sps[:, sps_index["actual_hy_Mg_Fe"]]*(1-oots[:, oots_index["Mg_num"]])*131.93) )#hy
    n_wt[:, nwt_index["mt"]] = np.where(condA,oots[:, oots_index["mt_Fe_III"]]*231.55   )                                          #mt
    n_wt[:, nwt_index["ilm"]] = np.where(condA,oots[:, oots_index["Ilm_Ti"]]*151.72   )                                           #ilm
    n_wt[:, nwt_index["ap"]] = np.where(condA,oots[:, oots_index["ap_P"]]*328.87     )                                         #ap
    nwt_columns = [index for col_name, index in nwt_index.items() if col_name != "total"]                       # Exclude "total"
    sum_of_nwt = np.where(condA,np.sum(n_wt[:, nwt_columns], axis=1))                                                                    # sum all columns except "residual_si"
    n_wt[:, nwt_index["total"]] = np.where(condA,sum_of_nwt)

    return n_wt, nwt_index

