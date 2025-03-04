###########################################################################################################################################################
#  Import libraries
###########################################################################################################################################################
from file_selector import FileSelector, load_data
import json

###########################################################################################################################################################
#  Intro messages
###########################################################################################################################################################
welcome_message = """
Welcome to epMin"""
instructions_message = """

This script reads a .DAT or .xlsx file that contains major and minor element contents of mineral phases.


Please bear in mind that the file should follow the output format from Probe Software:

1) Analyses need to be shown as rows
2) The heading of the file should have Sample id as "SAMPLE" and the element/oxides should not contain more than the chemical symbol 

"""
InstructionWindow = FileSelector(welcome_message, instructions_message) #creates window and gives the option to browse a file
file_path, file_type = InstructionWindow.get_file_path_and_type() #obtains the file type and path from the browse window

###########################################################################################################################################################
# check if file is valid and load the data
###########################################################################################################################################################
if file_type == "Unsupported file type.":
    print("The selected file is not supported.")
else:
    print(f"File selected: {file_path}")
    print(f"File type: {file_type}")

    # Load data using the appropriate class for the file type
    data = load_data(file_path, file_type)

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

class OxideDataCalc:
    def __init__(self, data, elements_weights, oxides_weights):
        self.data = data  # Data loaded from the file as a dictionary
        self.elements_weights = elements_weights  # Elements weight dictionary
        self.oxides_weights = oxides_weights  # Oxides weight dictionary

    def clean_data(self):
        clean_data = {}
        for sample_name, composition in self.data.items():
            total = composition.get('TOTAL', 0)  # Access TOTAL key instead of 'Total'
            feo = composition.get('FeO', 0)
            cao = composition.get('CaO', 0)

            # Flag data as bad if criteria are met
            if total < 80 and (feo < 40 or cao < 50):
                clean_data[sample_name] = {key: 0 for key in composition}
            else:
                clean_data[sample_name] = composition
        return clean_data

    def to_anhydrous(self):
        anhydrous_data = {}
        for sample_name, composition in self.data.items():
            total = composition.get('TOTAL', 1)  # Avoid division by zero

            # Normalizing to 100% if total is non-zero
            if total > 0:
                anhydrous_data[sample_name] = {
                    key: value * 100 / total if key != 'TOTAL' else 100 for key, value in composition.items()
                }
            else:
                anhydrous_data[sample_name] = composition
        return anhydrous_data

    def to_mol(self):
        molar_data = {}
        for sample_name, composition in self.data.items():
            total = composition.get('TOTAL', 1)  # Avoid division by zero
            molar_composition = {}

            # Convert composition to moles based on both element and oxide weights
            if total > 0:
                for key, value in composition.items():
                    if key in self.oxides_weights:
                        oxide_info = self.oxides_weights[key]  # Get oxide info from oxides_weights dictionary
                        molecular_weight = oxide_info['molecular_weight']
                        # Use molecular weight for mole conversion for oxides
                        molar_composition[key] = value / molecular_weight
                    elif key in self.elements_weights:
                        element_weight = self.elements_weights[key]  # Get element weight from elements_weights dictionary
                        # Use atomic weight for mole conversion for elements
                        molar_composition[key] = value / element_weight
                    else:
                        molar_composition[key] = value  # For anything else, just keep the value as is
                molar_composition['TOTAL'] = 100  # Ensure TOTAL is set to 100
            molar_data[sample_name] = molar_composition
        return molar_data

    def to_cat(self):
        cation_data = {}

        for sample_name, composition in self.data.items():
            total = composition.get('TOTAL', 1)  # Avoid division by zero
            cation_composition = {}

            # Calculate cation fractions for each oxide
            if total > 0:
                for key, value in composition.items():
                    if key != 'TOTAL':
                        # Check if key is an oxide and use the "non_oxygen_atoms" from the oxides_weights dictionary
                        if key in self.oxides_weights:
                            oxide_info = self.oxides_weights[key]
                            non_oxygen_atoms = oxide_info['non_oxygen_atoms']
                            cation_composition[key] = value * non_oxygen_atoms  # Apply non_oxygen_atoms for cation calculation

                        # Handle elements directly (without oxides) for cation composition
                        if key in self.elements_weights:
                            # Apply a simple logic, though element-specific cation fractions might need adjustments
                            cation_composition[key] = value  # Adjust according to element valence

                # Sum all cation fractions and store in 'TOTAL' column
                cation_composition['TOTAL'] = sum(cation_composition.values())

                # Normalize all cation fractions based on the sum
                for key in cation_composition:
                    if key != 'TOTAL':
                        cation_composition[key] /= cation_composition['TOTAL']
            cation_data[sample_name] = cation_composition
        return cation_data

    def add_fe2o3(self):
        Fe2_Fetot = 0.7  # Ratio Fe2+/Fe3+ for NNO
        updated_data = {}
        for sample_name, composition in self.data.items():
            feo = composition.get('FeO', 0)

            # Check if Fe2O3 already exists, if not, calculate and add it
            if 'Fe2O3' not in composition:
                Fe2O3 = feo * (1 - Fe2_Fetot) * 1.11134
                updated_composition = composition.copy()
                updated_composition['Fe2O3'] = Fe2O3
                updated_composition['FeO'] = feo * Fe2_Fetot
                updated_data[sample_name] = updated_composition
            else:
                updated_data[sample_name] = composition

        return updated_data
    
# 4. Initialize OxideDataCalc with the loaded data and weights
calc = OxideDataCalc(data, elements_weights, oxides_weights)

    # 5. Step 1: Add Fe2O3
data_clean = calc.clean_data()

print(data_clean)