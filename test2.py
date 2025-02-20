

import tkinter
from tkinter import ttk
from tkinter import filedialog

import sv_ttk

def get_file_path():
    global file_path
    # Open and return file path
    file_path= filedialog.askopenfilename(title = "Select A File", filetypes = (("text files", "*.DAT"), ("Excel spreadsheet", "*.xlsx")))
    l1 = ttk.Label(root, text = "File path: " + file_path).pack()
    sv_ttk.set_theme("dark")

welcome = """

Welcome to epMin
"""
message = """
epMin calculates mineral formulae and site activities from EPMA analyses.

Please be aware that this script is designed for the output file obtained with CalcImage (J. Donovan) when processing data. The input file for single-point analyses should be in .DAT or .xlsx format. The heading of the file should contain at least sample name as "SAMPLE" and the oxides as m/m%

Browse your file using the button below
"""

root = tkinter.Tk()
root.title("Select a file") #Heading of the window
root.geometry("600x350") #Window size
info = ttk.Label(root, text=welcome, wraplength=500, font=("Arial",12, "bold")).pack() #Label widget for the message with word wrapping
label = ttk.Label(root, text=message, wraplength=500, font=("Arial",12)).pack() #Label widget for the message with word wrapping
button = ttk.Button(root, text="Browse",command = get_file_path).pack()

# This is where the magic happens
sv_ttk.set_theme("dark")



root.mainloop()




