# %%
from joblib import Parallel, delayed
from tqdm import tqdm
import gc
import os
import shutil
# %%
folder = r"./scn_files/"
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
# %%


def Run(fl_info):
    fl, fln = fl_info
    os.system(f'python run.py --filename "{fl}" --path "{fln}"')


files_to_run = [(fl, f"{folder}/{fl}") for fl in get_only_files_w_xml]
Parallel(n_jobs=2)(delayed(Run)(itm) for itm in tqdm(files_to_run)) #You can increase the n_jobs based on your CPU resources
