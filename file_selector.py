###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
import tkinter as tk #for graphical user interface
from tkinter import filedialog, ttk #for browser and more
import sv_ttk #to make the window dark theme
import pandas as pd # to handle data coming as excel file
import csv #to handle data coming as a text file
import json
import os
import re

###########################################################################################################################################################
# to read files and set variables
###########################################################################################################################################################
class FileReader:  
    def __init__(self, file_path):  # constructor method
        self.file_path = file_path  # defines file_path with the input of the constructor
        self.data = {}  # initialises a dictionary to load the data in
    def read(self):  # empty public method to read
        raise NotImplementedError("forgot to implement method in subclass for .DAT or .xlsx")  # return error if somebody calls only this method
    def get_data(self):  # returns the data dictionary to the script that needs it
        return self.data

###########################################################################################################################################################
# subclass that will override the way the file is read (if tab-separated)
###########################################################################################################################################################
class DatFileReader(FileReader):  
    def read(self):  # public method to read .DAT file
        with open(self.file_path, newline='') as f: #opens the file using the path selected
            reader = csv.DictReader(f, delimiter='\t')  #Read as dictionary using tab as delimiter
            for row in reader: #reads line by line 
                sample_name = row['SAMPLE']  # "SAMPLE" is the column name for sample names
                self.data[sample_name] = extract_elementsoxides(row) # calls the function to only get compositional information
                self.data[sample_name]['TOTAL'] = float(row['TOTAL'])
###########################################################################################################################################################
# subclass that will override the way the file is read (if tab-separated)
###########################################################################################################################################################
class DatFileReaderMap(FileReader):  
    def read(self):  # public method to read .DAT file
        extraname = r'_FirstPass_\d{5}__Oxide_Image_Classify'  # Part of the file name to remove
        new_file_name = re.sub(extraname, '', os.path.basename(self.file_path))  # Simplify the file name
        directory = os.path.dirname(self.file_path)  # Directory of the input excluding the name
        new_file_path = os.path.join(directory, new_file_name)  # New path using the new name
        os.rename(self.file_path, new_file_path)  # Rename the original file
        print(f"File renamed to: {new_file_path}")
        with open(new_file_path, 'r') as file: #Fix the contents of the file
            lines = file.readlines()
        lines[1] = lines[1].replace('\t\t', '\t') # Correct the second line if there are two tab characters
        with open(new_file_path, 'w') as corrected_file:
            corrected_file.writelines(lines)
        with open(new_file_path, newline='') as f: # Now read the corrected file
            lines = f.readlines()
            header = lines[1].strip().split('\t')  # Use the second row as the header
            header = [col.replace('"', '').replace(" WT%", "").replace("wt%", "").strip() for col in header]
            data_lines = lines[2:]
            reader = csv.DictReader(data_lines, fieldnames=header, delimiter='\t')  # Read as dictionary using tab as delimiter
            for row in reader:  # Read line by line
                cleaned_row = {key.replace('"', '').strip(): value for key, value in row.items()}
                if 'X' in cleaned_row:
                    cleaned_row['X_coord'] = float(cleaned_row.pop('X'))
                if 'Y' in cleaned_row:
                    cleaned_row['Y_coord'] = float(cleaned_row.pop('Y'))

                NXY = cleaned_row['NXY']  # "NXY" is the column name for sample names
                self.data[NXY] = extract_elementsoxides(cleaned_row)  # Calls the function to only get compositional information
                self.data[NXY]['X_coord'] = float(cleaned_row['X_coord'])
                self.data[NXY]['Y_coord'] = float(cleaned_row['Y_coord'])
                self.data[NXY]['NX'] = int(cleaned_row['NX'])
                self.data[NXY]['NY'] = int(cleaned_row['NY'])
                self.data[NXY]['Total'] = float(cleaned_row['Total'])

###########################################################################################################################################################
# subclass that will override the way the file is read (Excel)
###########################################################################################################################################################
class XlsxFileReader(FileReader):
    def read(self): # public method to read
        df = pd.read_excel(self.file_path) #opens the file using the path selected
        for _, row in df.iterrows(): #reads line by line
            sample_name = row['SAMPLE']  # Assuming "SAMPLE" is the column name for sample names
            self.data[sample_name] = extract_elementsoxides(row) # calls the function to only get compositional information
            self.data[sample_name]['TOTAL'] = float(row['TOTAL'])

###########################################################################################################################################################
# subclass that will override the way the file is read (Excel)
###########################################################################################################################################################
class XlsxFileReaderMap(FileReader):
    def read(self): # public method to read
        df = pd.read_excel(self.file_path) #opens the file using the path selected
        for _, row in df.iterrows(): #reads line by line
            NXY = row['NXY']  # "SAMPLE" is the column name for sample names
            self.data[NXY] = extract_elementsoxides(row) # calls the function to only get compositional information
            self.data[NXY]['X'] = float(row['X'])
            self.data[NXY]['Y'] = float(row['Y'])
            self.data[NXY]['NX'] = int(row['NX'])
            self.data[NXY]['NY'] = int(row['NY'])

###########################################################################################################################################################
# main fuctions to load json files, to identify compositional columns in the data and to load the data
###########################################################################################################################################################
def load_json(): # function to load elements_weights and oxides_weights from JSON files and return two dictionaries
    with open('elements_weights.json', 'r') as f: #dictionary with elements and their atomic weights
        elements_weights = json.load(f)  
    with open('oxides_weights.json', 'r') as f: # dictionary with oxides and their molecular weights
        oxides_weights = json.load(f)
    return elements_weights, oxides_weights

elements_weights, oxides_weights = load_json()

def extract_elementsoxides(row):
        filtered_data = {} # empty dictionary
        for key in row.keys(): #loop over the columns in the data
            if key in elements_weights or key in oxides_weights: #compares columns to the dictionaries containing possible oxides or elements
                try:
                    filtered_data[key] = float(row[key])  # Convert to float
                except ValueError:
                    continue  # In case the value isn't a valid number
        return filtered_data

file_handlers = {
    'single-point': {
        '.DAT': DatFileReader,      # .DAT file for single-point data
        '.xlsx': XlsxFileReader    # .xlsx file for single-point data
    },
    'map': {
        '.DAT': DatFileReaderMap,   # .DAT file for map data
        '.xlsx': XlsxFileReaderMap  # .xlsx file for map data
    }
}


def load_data(file_path, file_type, data_structure): #fuction to load the data by selecting the appropriate method depending on file format
    try:
        reader_class = file_handlers[data_structure][file_type] #finds the type of reading according to the type
        reader = reader_class(file_path)
        reader.read()  # Read the data using the appropriate reader
        return reader.get_data() # gets data into dictionaries
    except KeyError:
        raise ValueError(f"Unsupported file type: {file_type}")
    
###########################################################################################################################################################
# main class that creates GUI and allows file browsing and file type check
###########################################################################################################################################################
class FileSelector:
    def __init__(self, titlew, message): #constructor method to select the file and load the data
        self.titlew = titlew # message title
        self.message = message # message to display to the user
        self.file_path = None #initialises the file path variable   
        self.file_type = None #initialises the file type variable
    def _create_window(self): #private method to design the window
        root = tk.Tk() # Create a Tkinter window
        root.title("Select a file") # Set window title
        myframe = tk.Frame(root)
        myframe.pack(fill='both',expand=True)
        ttk.Label(root, text=self.titlew, wraplength=500, font=("Arial", 12, "bold")).pack() # Display title
        ttk.Label(root, text=self.message, wraplength=500, font=("Arial", 12)).pack() # Display message
        button = ttk.Button(root, text="Browse", command=lambda: self._open_file_dialog(root)) # Browse button
        button.pack() # Display button
        sv_ttk.set_theme("dark") #set dark theme
        root.mainloop() #starts tkinter to do a loop
    def _open_file_dialog(self, root): #private method to keep the window but also open a browse window
        self.root = root #store reference for the window
        root.after(0, self._file_dialog) #open file browser
    def _file_dialog(self): # private method to design the browse window
        self.file_path = filedialog.askopenfilename(title="Select A File", filetypes=(("All files", "*.*"),))
        if self.file_path:
            self._validate_file_type() #calls the private method to check the file format
        self.root.destroy()  # close the window after file selection
    def _validate_file_type(self): #private method to check the type of the file that was loaded
        if self.file_path.endswith(".DAT"): #check if it is a text file
            self.file_type = ".DAT"
        elif self.file_path.endswith(".xlsx"): #check if it is excel file
            self.file_type = ".xlsx"
        else:
            self.file_type = "Unsupported file type." #raise an error
    def get_file_path_and_type(self): # public method to call the window and return the path selected and the type
        self._create_window()
        return self.file_path, self.file_type