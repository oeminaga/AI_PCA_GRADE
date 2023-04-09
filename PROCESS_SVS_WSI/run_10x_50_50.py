import argparse
from Objects import HistoImage
from plexusnet.architecture import LoadModel
import pandas as pd
import numpy as np
import tensorflow as tf
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
gpus = tf.config.list_physical_devices('GPU')

for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)
config = "10_50_50"
parser = argparse.ArgumentParser(description='get file')
parser.add_argument('--filename', help='filename')
parser.add_argument('--path', help='the full path')
args = parser.parse_args()

if __name__ == "__main__":
    fln = str(args.path)
    fl = str(args.filename)
    weight_file = "./weight_98.h5"
    model = LoadModel(weight_file)
    data = {"filename": [],
            "score": []}
    # REQUIRED MORE MEMORY AND MAY CAUSE OUT OF MEMORY ERROR - TESTED SUCCESSFULLY ON 128GB RAM
    X = HistoImage(fln, batch_size=16, target_magnification=10,
                   overlap_rate=(0.5, 0.5))
    P = model.predict(X, verbose=0, workers=4)
    fl_npy = fl.replace("svs", "npy")
    data["filename"].append(fl)
    data["score"].append(np.mean(P[:len(X.store_coordination)]))
    np.save(fr"../outcome/{config}/npy_prostate/{fl_npy}", P)
    fl_csv = fl.replace("svs", "csv")
    pd.DataFrame(data).to_csv(
        fr"../outcome/{config}/csv_prostate/{fl_csv}", index=False)
