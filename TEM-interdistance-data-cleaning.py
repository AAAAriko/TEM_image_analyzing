import pandas as pd
 
file_name = input("Name of the xlsx file, do not include .xlsx:")
threshold = int(input("Largest acceptable length:"))

sheet_fullname = file_name + ".xlsx"
cleaned_sheet_name = file_name + ".txt"

# read by default 1st sheet of an excel file
all_datasheets = pd.read_excel(sheet_fullname,sheet_name = None)
f = open(cleaned_sheet_name,'a')

for sheet in all_datasheets:
    lengths = all_datasheets[sheet].iloc[:,0]
    for length in lengths:
        if length <= threshold:
            valid_data = str(length)
            f.write(valid_data + "\n")

f.close()


