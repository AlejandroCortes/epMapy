import tkinter as tk

# Functions that may be activated based on the choices
def function_1(): print("Function 1 activated!")
def function_2(): print("Function 2 activated!")
def function_3(): print("Function 3 activated!")
def function_4(): print("Function 4 activated!")
def function_5(): print("Function 5 activated!")
def function_6(): print("Function 6 activated!")
def function_7(): print("Function 7 activated!")
def function_8(): print("Function 8 activated!")

options = [
    #("Plot K/Na", function_1),
    #("Plot NaK/Al", function_2),
    #("Plot NaKCa/Al", function_3),
    #("Plot zrc saturation", function_4),
    #("Plot normative mineralogy", function_5),
    ("Plot raw oxides",  function_5),
    #("Plot anhydrous oxides", function_7),
    ("save raw oxides as grids in json format",  function_7)
]
# A dictionary mapping options to functions and labels
def on_submit():
    selected_options = []
    for i, (label, func) in enumerate(options):
        if vars[i].get() == 1:
            selected_options.append(label)
            func()
    
    if selected_options:
        print(f"Running: {', '.join(selected_options)} ...")
    root.quit()  # Close the window

# Function to handle the form submission

# Set up the tkinter window
root = tk.Tk()
root.title("Select Options")

# Instruction label
instruction_label = tk.Label(root, text="Choose what you want to do. You can select multiple options:", font=("Arial", 12))
instruction_label.pack(pady=10)

# Create checkbuttons dynamically using a loop
vars = []
for label, func in options:
    var = tk.IntVar()
    vars.append(var)
    checkbutton = tk.Checkbutton(root, text=label, variable=var)
    checkbutton.pack(anchor='w')

# Submit button
submit_btn = tk.Button(root, text="Continue", command=on_submit)
submit_btn.pack(pady=10)

# Start the Tkinter event loop
root.mainloop()
