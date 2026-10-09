# Reference only: phylodeep's pretrained FFNN scalers (sklearn StandardScaler) as CSV rows: mean, then scale.
import sys, os, pickle, numpy as np
src, dst = sys.argv[1], sys.argv[2]
for f in sorted(os.listdir(src)):
    if f.endswith(".pkl"):
        sc = pickle.load(open(os.path.join(src, f), "rb"))
        np.savetxt(os.path.join(dst, f.replace(".pkl", ".csv")), np.vstack([sc.mean_, sc.scale_]), delimiter=",", fmt="%.17g")
        print(f, sc.mean_.shape)
