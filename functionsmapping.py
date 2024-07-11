#importing libraries for standalone functions
import numpy as np
import pandas as pd

def clean_data(X):
    a   = len(X[:,0]) #num of data points
    b   = len(X[0,:]) #num of oxides + total
    clean = np.copy(X)  #copies the input to avoid overwriting
    for i in range (a):
        if (clean[i,-1] <80) & ((clean[i,6] <40) | (clean[i,3] <50)):
            clean[i,:] = clean[i,:]*0 #normalising data to 100%
    return clean
#This function takes major oxide compositions and normalise them to 100%
def to_anhydrous(X):
    a   = len(X[:,0]) #num of data points
    b   = len(X[0,:]) #num of oxides + total
    anh = np.copy(X)  #copies the input to avoid overwriting
    for i in range (a):
        for j in range (b):
            if (anh[i,-1] > 0):
                anh[i,j] = anh[i,j]*100/anh[i,-1] #normalising data to 100%
    return anh
#This function takes major oxide compositions in wt.% and converts them to mol
def to_mol(X):
    a   = len(X[:,0]) #num of pixels
    b   = len(X[0,:]) #num of oxides + total
    #molecular weights with FeOt
    mw  = [60.08,79.87,101.96,71.85,70.94,40.3,56.08,61.98,94.2,141.94,100]
    #molecular weights with FeO and Fe2O3 calculated
    mw2 = [60.08,79.87,101.96,159.7,71.85,70.94,40.3,56.08,61.98,94.2,141.94,100]
    mol = np.copy(X)  #copies the input to avoid overwriting
    for i in range (a):
        for j in range (b):
            if mol[i,j]>0:
                if b == 11:
                    mol[i,j] = mol[i,j]/mw[j]  #converts to mol
                else:
                    mol[i,j] = mol[i,j]/mw2[j] #converts to mol
    return mol
#This function takes major oxide compositions in mol and converts them to cation fractions
def to_cat(X):
    a         = len(X[:,0])               #num of pixels
    b         = len(X[0,:])               #num of oxides + total
    ct        = [1,1,2,1,1,1,1,2,2,2,1]   #num of cations per oxide with FeOt
    ct2       = [1,1,2,2,1,1,1,1,2,2,2,1] #num of cations per oxide with FeO and Fe2O3 calc
    cat       = np.copy(X)                #copies the input to avoid overwriting
    for i in range (a):
        for j in range (b-1):
            if cat[i,-1]>0:
                if b == 11:
                    cat[i,j] = cat[i,j]*ct[j]  #calculate cation fractions
                else:
                    cat[i,j] = cat[i,j]*ct2[j] #calculate cation fractions
        if cat[i,-1]>0:
            cat[i,-1] = np.sum(cat[i,:-1]) #calculates the sum of oxide mol per analysis
            cat[i,:-1] = cat[i,:-1]/cat[i,-1]
    return cat

def add_fe2o3(X):
    Fe2_Fetot = 0.7         #ratio Fe2+/Fe3+Fe2+ for NNO
    X2        = np.copy(X)  #copies the input to avoid overwriting
    if len(X[0,:]) == 11:
        Fe2O3   = X[:,3]*(1-Fe2_Fetot)*1.11134 #calculate Fe2O3 from FeOtot and Fe2+/Fetot
        #add the Fe2O3 data to the array
        X2      = np.append(X2,np.reshape(Fe2O3,(len(Fe2O3),1)),axis=1)
        X2[:,3] = X2[:,3]*Fe2_Fetot #calculate FeO from FeOtot and Fe2+/Fetot ratio
    X3        = np.copy(X2) #copyng array to reshape
    X3[:,3]   = X2[:,11]    #locating Fe2O3 before FeO
    X3[:,4:]  = X2[:,3:-1]  #pasting the rest of the elements after Fe2O3
    return X3

def norm_calc(X):
    Z = np.copy(X)  #copies the input to avoid overwriting
    a = len(Z[:,0]) #num of pixels
    oots =np.zeros((a,23)) #calculating with oxides other than silica
    for i in range (a):
        oots[i,0]  = Z[i,10]*10/3                                       #ap Ca
        oots[i,1]  = Z[i,10]                                            #ap P
        oots[i,2]  = Z[i,1]                                             #ilm Fe(II)
        oots[i,3]  = Z[i,1]                                             #ilm Ti
        oots[i,4]  = Z[i,9]                                             #or/lc Al
        oots[i,5]  = Z[i,9]                                             #or/lc K
        if Z[i,7] > (Z[i,2]-Z[i,5]):
            oots[i,6] = Z[i,2]-Z[i,5]               #ab/ne Al
        else:
            oots[i,6] = Z[i,8]                      #ab/ne Al
        oots[i,7]  = oots[i,6]                                          #ab/ne Na
        if (Z[i,7]-oots[i,7]) < Z[i,3]:
            oots[i,9] = Z[i,8]-oots[i,7]            #ac Fe(III)
        else:
            oots[i,9]=Z[i,3]                        #ac Fe(III)
        oots[i,8]  = oots[i,9]                                          #ac Na
        oots[i,10] = Z[i,8]-oots[i,7]-oots[i,8]                         #ns Na
        oots[i,12] = Z[i,3]-oots[i,9]                                   #mt Fe(III)
        oots[i,11] = oots[i,12]                                         #mt Fe(II)
        if oots[i,8] > 0.000001:
            oots[i,13] = 0                          #an Al prelim
        else:
            oots[i,13] = Z[i,2]-oots[i,4]-oots[i,6] #an Al prelim
        if oots[i,13] > (Z[i,7]-oots[i,0]):
            oots[i,14] = Z[i,7]-oots[i,0]           #an Al final
        else:
            oots[i,14] = oots[i,13]                 #an Al final
        if oots[i,13] == oots[i,14]:
            oots[i,15] = 0                          #C Al
        else:
            oots[i,15] = oots[i,13]-oots[i,14]      #C Al
        oots[i,16] = oots[i,14]                                         #an Ca
        oots[i,17] = Z[i,7]-oots[i,0]-oots[i,16]                        #di Ca
        if Z[i,-1]>0:
            oots[i,22] = Z[i,6]/(Z[i,6]+Z[i,4]+Z[i,5]-oots[i,2]-oots[i,11]) #XMg
        oots[i,18] = oots[i,17]*oots[i,22]                              #di Mg
        oots[i,19] = oots[i,17]*(1-oots[i,22])                          #di FeII 
        oots[i,20] = Z[i,4]+Z[i,5]-oots[i,2]-oots[i,11]-oots[i,19]      #ol/hy Fe(II)
        oots[i,21] = Z[i,6]-oots[i,18]                                  #ol/hy Mg

    fps=np.zeros((a,8)) #First pass silica
    for i in range (a):
        fps[i,0]   = 4*oots[i,5]                                        #lc
        fps[i,1]   = 2*oots[i,7]                                        #ne
        fps[i,2]   = 4*oots[i,8]                                        #ac
        fps[i,3]   = oots[i,10]                                         #ns
        fps[i,4]   = 2*oots[i,16]                                       #an
        fps[i,5]   = 2*oots[i,17]                                       #di
        fps[i,6]   = (oots[i,20]+oots[i,21])/2                          #ol
        #residual Si
        fps[i,7]   = Z[i,0]-(fps[i,0]+fps[i,1]+fps[i,2]+fps[i,3]+fps[i,4]+fps[i,5]+fps[i,6])
  
    sps=np.zeros((a,11)) #Second pass silica
    for i in range (a):
        if (oots[i,5]*2) > fps[i,7]:
            sps[i,0]  = (fps[i,0]+(fps[i,7]-(oots[i,5]*2)))/2 #or K
        else:
            sps[i,0]  = oots[i,5]                             #or K
        sps[i,1] = oots[i,5]-sps[i,0]                         #lc K
        sps[i,2] = fps[i,7]-(2*sps[i,0])                      #residual Si
        if (oots[i,7]*4) > sps[i,2]:
            sps[i,3]  = (sps[i,2]+fps[i,1]-(2*oots[i,7]))/4   #ab Na
        else:
            sps[i,3]  = oots[i,7]                             #ab Na
        sps[i,4] = oots[i,7]-sps[i,3]                         #ne Na
        sps[i,5] = sps[i,2]-(4*sps[i,3])                      #residual Si
        if sps[i,4] > 0.0000001:
            sps[i,6]  = oots[i,20]+oots[i,21]                 #prelim ol Mg+Fe
        else:
            sps[i,6]  = oots[i,20]+oots[i,21]-(2*sps[i,5])    #prelim ol Mg+Fe
        sps[i,7] = oots[i,20]+oots[i,21]-sps[i,6]             #prelim hy Mg+Fe
        if sps[i,6] < 0:
            sps[i,8]  = 0                                     #actual ol Mg+Fe
            sps[i,9]  = oots[i,20]+oots[i,21]                 #actual hy Mg+Fe
            sps[i,10] = sps[i,5]-((oots[i,20]+oots[i,21])/2)  #residual Si (Q)
        else:
            sps[i,8]  = sps[i,6]                              #actual ol Mg+Fe
            sps[i,9]  = sps[i,7]                              #actual hy Mg+Fe
            sps[i,10] = 0                                     #residual Si (Q)

    fs=np.zeros((a,7)) #final silica
    for i in range (a):
        fs[i,0]  = 6*sps[i,0]   #or
        fs[i,1]  = 4*sps[i,1]   #lc
        fs[i,2]  = 6*sps[i,3]   #ab
        fs[i,3]  = 2*sps[i,4]   #ne
        fs[i,4]  = 0.5*sps[i,8] #ol
        fs[i,5]  = sps[i,9]     #hy
        fs[i,6]  = sps[i,10]    #Q
    
    n_wt=np.zeros((a,16)) #recalculate to normative mineral in wt%
    for i in range(a):
        n_wt[i,0]  = fs[i,6]*60.08                                                 #Q
        n_wt[i,1]  = sps[i,0]*556.64                                               #or
        n_wt[i,2]  = sps[i,1]*436.48                                               #lc
        n_wt[i,3]  = sps[i,3]*524.42                                               #ab
        n_wt[i,4]  = sps[i,4]*284.1                                                #ne
        n_wt[i,5]  = oots[i,16]*278.2                                              #an
        n_wt[i,6]  = oots[i,15]*101.96                                             #C
        n_wt[i,7]  = oots[i,8]*462                                                 #ac
        n_wt[i,8]  = oots[i,10]*122.06                                             #ns
        n_wt[i,9]  = (oots[i,17]*116.16)+(oots[i,18]*100.38)+(oots[i,19]*131.39)   #di
        n_wt[i,10] = (sps[i,8]*70.34*oots[i,22])+(sps[i,8]*(1-oots[i,22])*101.89)  #ol
        n_wt[i,11] = (sps[i,9]*100.38*oots[i,22])+(sps[i,9]*(1-oots[i,22])*131.93) #hy
        n_wt[i,12] = oots[i,12]*231.55                                             #mt
        n_wt[i,13] = oots[i,3]*151.72                                              #ilm
        n_wt[i,14] = oots[i,1]*328.87                                              #ap
        n_wt[i,15] = np.sum(n_wt[i,:-1])                                           #total
    return n_wt

