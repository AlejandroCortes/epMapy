import tkinter as tk
import json
import sv_ttk

# Read JSON data from a file
with open('minerals.json', 'r') as f:
    minerals_data = json.load(f)  # This loads the content of the file directly

def on_select():
    print(f"Selected: {selected_value.get()}")

# Create the main window
root = tk.Tk()
root.title("Minerals Selection")

# Apply dark theme using sv_ttk
sv_ttk.set_theme("dark")

# Set a larger window size to fit the content better
root.geometry("500x450")  # Increased window size

# Add some styling (font and padding)
root.config(bg="#2e2e2e")  # Dark background for the window

# Create a Tkinter variable to hold the selected radio button value
selected_value = tk.StringVar(value="Olivine")  # Set the default value to 'Olivine' so it appears selected

# Create a frame to hold the radio buttons for better performance and layout
frame = tk.Frame(root, bg="#2e2e2e")  # Dark frame background
frame.pack(pady=20, padx=20, fill="both", expand=True)

# Add the instruction message above the radio buttons
instruction_label = tk.Label(frame, text="Please select a mineral from the list below:", font=("Arial", 14), fg="#ffffff", bg="#2e2e2e")
instruction_label.grid(row=0, column=0, columnspan=2, pady=10)

# Create two columns of radio buttons
columns = 2  # Set the number of columns
row = 1  # Start the radio buttons below the instruction label
column = 0

# Loop to create radio buttons in two columns for minerals from the JSON file
for index, mineral in enumerate(minerals_data.keys()):
    radio_button = tk.Radiobutton(frame, text=mineral, value=mineral, variable=selected_value, command=on_select, font=("Arial", 12), fg="#ffffff", bg="#2e2e2e", selectcolor="#4a4a4a")
    
    # Place the radio buttons in two columns
    radio_button.grid(row=row, column=column, sticky="w", pady=5, padx=10)
    
    # Update the row and column positions for the next radio button
    column += 1
    if column == columns:
        column = 0
        row += 1

# Add the "Multiple" option as an additional radio button
multiple_radio_button = tk.Radiobutton(frame, text="Multiple", value="Multiple", variable=selected_value, command=on_select, font=("Arial", 12), fg="#ffffff", bg="#2e2e2e", selectcolor="#4a4a4a")
multiple_radio_button.grid(row=row, column=column, sticky="w", pady=5, padx=10)

# Run the application
root.mainloop()