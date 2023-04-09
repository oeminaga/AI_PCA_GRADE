# %%
import pandas as pd
import gc
import os
import shutil
# %%
folder = r"../../Diagnostic_image/"

files_list = {}
for dirpath, dirnames, filenames in os.walk(folder):
    check_files_xml = []
    check_files_svs = []
    for filename in filenames:
        if filename.endswith(".xml"):
            check_files_xml.append(filename)
        if filename.endswith(".svs"):
            check_files_svs.append(filename)
    if len(check_files_xml) == 1 and len(check_files_svs) == 1:
        files_list[f"{dirpath}/{check_files_svs[0]}"] = f"{dirpath}/{check_files_xml[0]}"
# %%
len(files_list)
# %%
ID_file = []
ID_Replicate = []
ID_FILENAME = []
CASE_ID = []
Filename_path = []
Filename_xml = []
for key in files_list:
    filename = os.path.basename(key)
    ID_s = filename.split(".")[0]
    ID_file.append(ID_s)
    ID_FILENAME.append(filename)
    Filename_path.append(key)
    Filename_xml.append(files_list[key])
pd.DataFrame({
    "CASE_ID": ID_file,
    "FilenameID": ID_FILENAME,
    "Filename": Filename_path,
    "Annotation": Filename_xml

}).to_csv("svs_files.csv", index=False)

# %%
