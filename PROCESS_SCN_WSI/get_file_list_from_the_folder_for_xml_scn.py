# %%
import pandas as pd
import gc
import os
import shutil
# %%
folder = r"../HE_images/"
files_svs = [f for f in os.listdir(folder) if f.endswith("scn")]
# %%
get_only_files_w_xml = []
for fl in files_svs:
    fl_xml = fl.replace("scn", "xml")
    if os.path.exists(f"{folder}/{fl_xml}"):
        get_only_files_w_xml.append(fl)

# %%
print("len(get_only_files_w_xml)", len(get_only_files_w_xml))
# %%
ID_file = []
ID_Replicate = []
ID_FILENAME = []
CASE_ID = []
for itm in get_only_files_w_xml:
    ID_s = itm.split(".")[0]
    info_s = ID_s.split("-")
    rep = 0
    if len(info_s) > 1:
        rep = int(info_s[1])
    ID_file.append(ID_s)
    ID_Replicate.append(rep)
    CASE_ID.append(info_s[0])
    ID_FILENAME.append(itm)
pd.DataFrame({
    "Filename": ID_FILENAME,
    "FILENAME_ID": ID_file,
    "CASE_ID": CASE_ID,
    "Replicate": ID_Replicate

}).to_csv("scn_files.csv", index=False)

# %%
