# %%
import pandas as pd
import os
# %%
outcomes = os.listdir("../outcome")

for outcome in outcomes:
    source_ = f"../outcome/{outcome}/csv_prostate/"
    data_collection = []
    for file in os.listdir(source_):
        data = pd.read_csv(f"{source_}/{file}")
        data_collection.append(data)
    database = pd.concat(data_collection)
    itms = database.filename.str.split(".", expand=True)
    itms = itms[0].str[:15]  # .split("-", n=4,expand=True)
    database["ID_SAMPLE"] = itms
    database[outcome] = database["score"]
    database.to_csv(f"./{outcome}.csv", index=False)
