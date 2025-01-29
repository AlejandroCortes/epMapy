import numpy as np #math functions as pi,e,random...
import matplotlib.pyplot as plt #default plotting library
import pandas as pd #for printing nice tables (pd.DataFrame)
import seaborn as sbn #extra version of plt
import matplotlib.patches as mpatches #from brokenaxes import brokenaxes
import matplotlib
from matplotlib.path import Path
from matplotlib.transforms import Bbox

# Reading excel file with major and trace elements in bulk rocks
data = pd.read_excel(r'C:\\Users\\alejc\\Documents\\GitHub\\epMapy\\Finley.xlsx')

# List all oxide columns (excluding 'Instrument' and 'Phase' columns)
oxide_columns = [col for col in data.columns if col not in ['Instrument', 'Phase', 'Type', 'SAMPLE', 'TOTAL']]

# Loop through each sample
samples = data['SAMPLE'].unique()

for sample in samples:
    sample_data = data[data['SAMPLE'] == sample]
    print(sample_data)

    # Set up the figure for subplots (3 columns, enough rows for all oxides)
    num_oxides = len(oxide_columns)
    num_cols = 3  # Number of columns
    num_rows = (num_oxides + num_cols - 1) // num_cols  # Calculate rows needed

    # Create subplots with a specific number of rows and columns
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 5))

    # Flatten axes array for easier indexing
    axes = axes.flatten()
    errors = [0.00156787,0.00972674,0.00313189,0.00733235,0.107725,0.00704432,0.00464168,0.0101844,0.00916893,0.0650152,0.744531,0.0999656,0.365345]

    # Loop over the oxides and plot on the corresponding subplot
    for i, oxide in enumerate(oxide_columns):
        concentrations = sample_data[oxide]
        instruments = sample_data['Instrument']
        instruments2 = instruments[instruments !='Reference']
        if 'secondary' in sample_data['Type'].values:
            lim_1 = sample_data[sample_data['Instrument'] == "Reference"][oxide] + 1.5
            lim_2 = sample_data[sample_data['Instrument'] == "Reference"][oxide] - 1.5
            axes[i].scatter(instruments, concentrations, color='red')
        else:
            # If 'Type' is not 'secondary', you can just set arbitrary limits or leave it
            lim_1 = concentrations.max() + 0.1  # Just using the max concentration if not secondary
            lim_2 = concentrations.min() - 0.1  # Just using the min concentration if not secondary
            axes[i].scatter(instruments2, concentrations, color='blue')

        # Scatter plot for each oxide on the corresponding subplot
        
        axes[i].set_title(f'{oxide} ({sample})')
        axes[i].set_ylim(lim_2.min(), lim_1.max())  # Set the y-limits based on the reference or max/min
        axes[i].set_ylabel(f'{oxide} m/m%')
        axes[i].set_xlabel('Instrument')

    # Remove any unused subplots if oxides are fewer than the number of subplots
    for j in range(num_oxides, len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure to a file (change the filename as needed)
    plt.savefig(f'{sample}_oxides_plot.pdf', dpi=300)
    plt.close()  # Close the figure to prevent it from displaying immediately