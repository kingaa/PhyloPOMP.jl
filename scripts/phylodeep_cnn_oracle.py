# Reference only (Python is never the pipeline path): raw outputs of a phylodeep pretrained CNN (e.g. BD_SMALL_CNN.h5)
# on the rows of a CSV from phylodeep_cblv_oracle.py.
# usage: venv/bin/python -I scripts/phylodeep_cnn_oracle.py <model.h5> <cblv.csv> <out.csv>
import sys, os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import numpy as np

from tensorflow import keras
model = keras.models.load_model(sys.argv[1], compile=False)
X = np.loadtxt(sys.argv[2], delimiter=",")[:, :-1]
np.savetxt(sys.argv[3], model.predict(X, verbose=0).astype(float), delimiter=",", fmt="%.10g")
print(X.shape, "->", sys.argv[3])
