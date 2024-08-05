import os
import datetime
year = f'{datetime.date.year}+/'
instr = input('Are you using the JEOL 8530F? (Y/N):     ')
if instr == 'Y':
    instr = 'JEOL_JXA_8530F_Plus/'
    print('Setting up for JEOL JXA 8530F Plus ...       ')
else:
    instr = 'Cameca_SX100_EPMA'
    print('Setting up for Cameca SX100 ...      ')

name = input('Please type the name and surname of the user as follows Name_Surname:     ')
if os.path.exists('//fls-sip/share/workspace/iac-workspace/Probes/'):
    print("the path exists")
else:
    print("foolish")
print ("hello")
print ("hello")
print ("hello")
print ("hello")
print ("hello")
