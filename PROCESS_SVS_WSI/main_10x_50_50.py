# %% Load libraries
import gc
import os
import shutil
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd
# %% Load the configuration file, create the output and the list of svs files
config = "10x_50_50"
# This is the folder that contains the csv files for patchwise analysis per slide
os.makedirs(f"../outcome/{config}/csv_prostate", exist_ok=True)
# This is the folder that contains the npy files for patchwise analysis per slide
os.makedirs(f"../outcome/{config}/npy_prostate", exist_ok=True)
# This is the csv file that contains the list of svs files and corresponds to the job list
data = pd.read_csv("svs_files.csv")
# %%


def Run(fl_info):
    fl, fln = fl_info
    os.system(f'python run_{config}.py --filename "{fl}" --path "{fln}"')


# %% Apply parallel computing to run the code
files_to_run = [[x, y] for (x, y) in zip(
    data["CASE_ID"].to_list(), data["Filename"].to_list())]
Parallel(n_jobs=6)(delayed(Run)(itm) for itm in tqdm(files_to_run))
