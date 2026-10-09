# Reference only. Shows (1) that the unattached `keras.layers.Dropout(0.5)` pattern of the deeptimelearning repository
# (Lambert, Voznica & Morlon, diversification paper; not the 2022 PhyloDeep training code) builds no Dropout layer, and
# (2) that phylodeep's released BD_SMALL_CNN contains no dropout. It does not show how the 2022 PhyloDeep models were trained.
# usage: venv/bin/python -I scripts/phylodeep_dropout_check.py <phylodeep pretrained_models/models dir> <scratch dir>
import os, sys, json
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"; os.environ["CUDA_VISIBLE_DEVICES"] = ""
import numpy as np, h5py, tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Dense

# 1. Voznica's pattern (BD_ffnn_SS_mae.py:178-186): bare keras.layers.Dropout(0.5), no model.add
m = keras.Sequential()
m.add(keras.Input(shape=(98,)))
m.add(Dense(64, activation='elu')); keras.layers.Dropout(0.5)
m.add(Dense(32, activation='elu')); keras.layers.Dropout(0.5)
m.add(Dense(16, activation='elu')); keras.layers.Dropout(0.5)
m.add(Dense(8, activation='elu'));  keras.layers.Dropout(0.5)
m.add(Dense(4, activation='linear'))
print("1. Voznica pattern, layers:", [type(l).__name__ for l in m.layers])

# 2. Same model with model.add(Dropout): Dropout is saved in the .h5 config
m2 = keras.Sequential([keras.Input(shape=(98,)), Dense(64, activation='elu'), keras.layers.Dropout(0.5), Dense(4)])
m2.save(sys.argv[2] + "/with_dropout.h5")
with h5py.File(sys.argv[2] + "/with_dropout.h5") as f:
    cfg = json.loads(f.attrs["model_config"])
print("2. attached Dropout, saved .h5 layers:", [l["class_name"] for l in cfg["config"]["layers"]])

# 3. phylodeep's released BD_SMALL_CNN: training=True (dropout active if present) gives identical outputs
model = keras.models.load_model(sys.argv[1] + "/BD_SMALL_CNN.h5", compile=False)
x = np.random.default_rng(0).random((5, 402)).astype("float32")
outs = [model(x, training=True).numpy() for _ in range(20)]
print("3. released CNN, 20 calls with training=True, max spread:", max(np.abs(o - outs[0]).max() for o in outs))
m2out = [m2(x[:, :98], training=True).numpy() for _ in range(20)]
print("   control (model WITH dropout), same test, max spread:", max(np.abs(o - m2out[0]).max() for o in m2out))
