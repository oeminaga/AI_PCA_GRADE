import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="1"
import tensorflow as tf
import numpy as np
import pandas as pd
gpus = tf.config.list_physical_devices('GPU')

for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)
from plexusnet.architecture import LoadModel
from Objects import HistoImage
import argparse

parser = argparse.ArgumentParser(description='get file')
parser.add_argument('--filename', help='filename')
parser.add_argument('--path', help='the full path')

args = parser.parse_args()

if __name__=="__main__":
    fln=str(args.path)
    fl=str(args.filename)
    weight_file="./weight_98.h5"
    model=LoadModel(weight_file)
    data={"filename":[],
    "score":[]}
    X=HistoImage(fln, batch_size=16, overlap_rate=(0.5,0.5),target_magnification=10)
    P=model.predict(X, verbose=0, workers=4)
    fl_npy=fl.replace("scn", "npy")
    data["filename"].append(fl)
    data["score"].append(np.mean(P[:len(X.store_coordination)]))
    np.save(fr"../npy_prostate/{fl_npy}", P)
    fl_csv=fl.replace("scn", "csv")
    pd.DataFrame(data).to_csv(fr"../csv_prostate/{fl_csv}", index=False)