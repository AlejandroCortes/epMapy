import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import pickle
from scipy.spatial import distance

# Load the data
with open('quant_maps.pickle', 'rb') as file:
    data = pickle.load(file)

# Extract Li concentration map
li_conc_array = np.array(data['Li ([7Li]+)_conc'])

# Set colormap with vmax = 40
cmap = cm.get_cmap('viridis').copy()
cmap.set_under('black')
cmap.set_over('black')
norm = colors.Normalize(vmin=0, vmax=40, clip=False)

# Show image and get two clicked points
plt.imshow(li_conc_array, cmap=cmap, norm=norm)
plt.title('Click two points to extract a profile')
plt.axis('on')
points = plt.ginput(2, timeout=-1)
plt.close()

# Extract clicked coordinates
start = np.array(points[0])
end = np.array(points[1])
PIXEL_SIZE_UM = 5  # microns per pixel
STEP_SIZE_UM = 15  # we want samples every 10 microns
STEP_PIXELS = STEP_SIZE_UM / PIXEL_SIZE_UM  # = 2 pixels

line_length_px = distance.euclidean(start, end)
line_length_um = line_length_px * PIXEL_SIZE_UM  # now in microns

num_steps = int(line_length_um // STEP_SIZE_UM)

# Then interpolate that many points
x_vals = np.linspace(start[0], end[0], num_steps)
y_vals = np.linspace(start[1], end[1], num_steps)

# Extract profile values
profile = []
for x, y in zip(x_vals, y_vals):
    x_int, y_int = int(round(x)), int(round(y))
        # Define neighborhood bounds
    y_min = max(0, y_int - 1)
    y_max = min(li_conc_array.shape[0], y_int + 2)
    x_min = max(0, x_int - 1)
    x_max = min(li_conc_array.shape[1], x_int + 2)

    # Extract the 3x3 neighborhood
    neighborhood = li_conc_array[y_min:y_max, x_min:x_max]
    
    # Compute the mean, ignoring NaNs if any
    mean_value = np.nanmean(neighborhood)
    profile.append(mean_value)
profile = np.array(profile)
distance_um = np.arange(num_steps) * STEP_SIZE_UM

# === Improved Scatter Plot of the Profile ===
fig, ax = plt.subplots(figsize=(7, 6))  # Wider aspect ratio for better readability
ax.scatter(distance_um, profile, color='blue', s=60)  # Increased marker size
ax.errorbar(distance_um, profile, yerr= profile*0.1, mfc='blue', ms= 60, linestyle='')

# Increase font sizes
ax.set_xlabel(r'Distance ($\mu$m)', fontsize=14)
ax.set_ylabel(r'Li Concentration ($\mu$g/g)', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.grid(False)

plt.tight_layout()
plt.savefig('li_profile_scatter_10um.png', dpi=600)  # Higher resolution (600 DPI)
print("Profile scatter plot saved as 'li_profile_scatter_10um.png'.")

# Option to save the colormap with vmax=40
save_option = input("Do you want to save the colormap with vmax=40? (y/n): ").lower()
if save_option == 'y':
    plt.figure()
    plt.imshow(li_conc_array, cmap=cmap, norm=norm)
    plt.colorbar(label=r'Li ($\mu$g/g)')
    plt.axis('off')
    plt.savefig('li_concentration_cmap_vmax40.png', dpi=300, bbox_inches='tight')
    print("Colormap saved as 'li_concentration_cmap_vmax40.png'.")
