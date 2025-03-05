# epMapy

This repository provides tools for plotting EPMA quantitative maps and performing calculations on the data. It is specifically designed for processing output files obtained from CalcImage (J. Donovan).

The input file for this script must be in `.xlsx` format. Therefore, if your data is in `.DAT` format, you will need to convert it to `.xlsx` before using this script.

If you are using software other than CalcImage to generate the maps, please ensure that the file's header contains at least the following column names (case-sensitive):

X, Y: EPMA position coordinates (x-y).
NX, NY: Pixel number in the x and y directions (matrix dimensions).
NXY: Consecutive pixel number.
- The remaining columns contain weight percentages, including the Total column.

Example columns in the file:

`X, Y, NX, NY, NXY, SiO2 WT%, TiO2 WT%, Al2O3 WT%, FeO WT%, MnO WT%, MgO WT%, CaO WT%, Na2O WT%, K2O WT%, P2O5 WT%, Total`

## Getting Started

### Step 1: Create a Virtual Environment

You only need to do this once per machine:

`python -m venv .venv`

### Step 2: Activate the Virtual Environment

Always activate the virtual environment before running `python`, `pip`, etc. The prompt should begin with `(.venv)`:

For Windows (PowerShell):

```
Set-ExecutionPolicy Unrestricted -Scope Process
.\.venv\Scripts\activate
```

## Main Files

`epMap`: Used to call and run individual functions to plot maps from the probe data.
`epMin`: Used to call and run individual functions for recalculating mineral formulae from single-point analyses.

## Installation

To install the required dependencies, run the following commands in your terminal:

tkinter
```
pip install tkinter
```

sv_ttk
```
pip install sv_ttk
```

pandas
```
pip install pandas
```

numpy
```
pip install numpy
```

matplotlib
```
pip install matplotlib
```

scipy.ndimage
```
pip install scipy.ndimage
```

openpyxl
```
pip install openpyxl
```

## Function Descriptions

`show_intro()`: Displays instructions and general information about the script.

`browse_win()`: Opens a window to browse for files.

`clean_data()`: Cleans the data by removing low totals (e.g., epoxy, vesicles, cracks).

`to_anhydrous()`: Calculates the anhydrous base.

`to_mol()`: Converts weight percentages to molar proportions.

`to_cat()`: Converts molar to cation proportions.

`add_fe2o3()`: Calculates Fe2O3 from FeO total for normative mineralogy.

`norm_calc()`: normative mineralogy calculations.

`extract_profile()`: Extracts compositional profiles from the maps.

`load_json()`: Loads oxide and element weights, charges, etc.

`read()`: Reads `.txt` or `.xlsx` files containing compositional information.

`get_data()`: Returns the data dictionary to the calling script.

`extract_elementsoxides()`: Compares columns to dictionaries containing possible oxides or elements.

`load_data()`: Builds dictionaries; currently implemented only for the `epMin` script.

`_create_window()`: A private method to design the window.

`_open_file_dialog()`: A private method to open a browse window while retaining the current window.

`_file_dialog()`: A private method for designing the browse window.

`_validate_file_type()`: A private method to check the file type of the loaded file.

`get_file_path_and_type()`: A public method to open the browse window and return the selected file path and type.