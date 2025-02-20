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
import periodictable

###########################################################################################################################################################
# Important information for User
###########################################################################################################################################################
Introd = """
epMin calculates mineral formulae and site activities EPMA quantitative maps and performing calculations with the data within, 
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

print(periodictable.H.mass)
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

