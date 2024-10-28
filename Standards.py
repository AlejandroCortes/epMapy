import pandas as pd
import os
import sys

# Database Name
file_name = 'EPMA_Standards.csv'
# Database basic arguments
basic_arguments = ['Name', 'ID_Old', 'BM', 'Type', 'Block', 'Locality', 'Reference']
# Check if the CSV file already exists
if os.path.exists(file_name):
    # Load existing data if file exists
    data = pd.read_csv(file_name)
else:
    # Create a new DataFrame with the initial columns if file doesn't exist
    data = pd.DataFrame(columns=basic_arguments)

def add_entry():
    global data
    # Dictionary to store entry values
    entry = {}
    # Entry values for initial columns only
    for col in basic_arguments:
        entry[col] = input(f"Enter value for {col}: ")
    # Ask for element composition
    while True:
        add_argument = input("Do you want to add an element? (y/n): ")
        if add_argument.lower() == 'y':
            # Ask which element
            new_argument = input("Enter the name of the element: ")
            new_value = input(f"Enter value for {new_argument}: ")

            # Add the new argument to the entry and DataFrame if it doesn't exist
            entry[new_argument] = new_value
            if new_argument not in data.columns:
                data[new_argument] = pd.NA  # Add a new column with empty values in existing entries
        else:
            break

    # Convert the entry dictionary to a DataFrame and concatenate with the existing DataFrame
    new_entry_df = pd.DataFrame([entry])
    data = pd.concat([data, new_entry_df], ignore_index=True)
    print("Entry added successfully!")

def modify_entry():
    global data

    # Display existing entries
    print("\nCurrent Entries:")
    print(data)

    # Prompt for the index of the entry to modify
    index = int(input("Enter the index of the entry you want to modify (0 to {}): ".format(len(data) - 1)))

    if index < 0 or index >= len(data):
        print("Invalid index. No modifications made.")
        return

    # Display the selected entry
    print("\nSelected Entry:")
    print(data.iloc[index])

    # Prompt for which field to modify
    field_to_modify = input("Enter the name of the field you want to modify: ")

    if field_to_modify in data.columns:
        new_value = input(f"Enter the new value for {field_to_modify}: ")
        data.at[index, field_to_modify] = new_value
        print("Entry updated successfully!")
    else:
        print("Invalid field name. No modifications made.")

def main():
    print("NHM EPMA Standard Database")
    
    while True:
        
        action = input("\nChoose an action: (1) Add Entry (2) Modify Entry (3) Quit: ")

        if action == '1':
            add_entry()
        elif action == '2':
            modify_entry()
        elif action == '3':
            break
        else:
            print("Invalid option, please try again.")

    # Sort the DataFrame by 'field3' alphabetically before saving
    data.sort_values(by='Type', inplace=True, ignore_index=True)

    # Save DataFrame to CSV
    data.to_csv(file_name, index=False)
    print(f"All entries have been saved to {file_name}, sorted alphabetically by 'Type'.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted. Saving data...")
        
        # Sort by 'Type' alphabetically before saving
        data.sort_values(by='Type', inplace=True, ignore_index=True)
        
        data.to_csv(file_name, index=False)
        print(f"Data saved to {file_name}, sorted alphabetically by 'Type'.")
        sys.exit()

