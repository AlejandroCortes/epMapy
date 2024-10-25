import pandas as pd
import glob

file_path = r'\\fls-sip\share\workspace\iac-workspace\Probes\JEOL_JXA_8530F_Plus\2024\Beth_Jones\241023_Jones\Glass'

def compile_files(file_path):
    txt_files = glob.glob(file_path + '/*.csv')
    df_list = (pd.read_csv(file, sep = ',') for file in txt_files)
    compiled_df = pd.concat(df_list, ignore_index = True)
    compiled_df.to_csv(file_path + '\\241025_Jones.csv')
    #compiled_df.columns.values[0] = 'Sample'
    return compiled_df

print (compile_files(file_path))