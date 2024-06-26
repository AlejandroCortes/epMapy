epMapy.py
===============

This repository allows plotting EPMA quantitative maps and do calculations with the data within.

Please be aware that this script is designed for the output file obtained with CalcImage (J. Donovan).
The input file for this script is a .xslx file, so you should convert .DAT to .xslx before

If another software is used to obtain the maps make sure that the file's head contains 
at least the following column names (case and character sensitive):
X,Y,NX,NY,NXY,SiO2 WT%,TiO2 WT%,Al2O3 WT%,FeO WT%,MnO WT%,MgO WT%,CaO WT%,Na2O WT%,K2O WT%,P2O5 WT%,Total
-X and Y: x-y coordinates EPMA position
-NX and NY: pixel number in x-y (matrix)
-NXY: consecutive pixel number, all the rest values are weight percents including the total. 

Unfortunatelly the scripts do not work with more oxides or halogens at the moment, but can be easily
modified to do so by the user. Just keep in mind to also modify the external functions that are being used.