# Reference only (Python is never the pipeline path): phylodeep's encode_into_summary_statistics per tree (99 values +
# rescale factor), and, when the pretrained_models dir is given, BD_SMALL_FFNN's scaler (mean, scale) and raw outputs.
# usage: venv/bin/python -I scripts/phylodeep_ss_oracle.py <in.tsv: p_sample, plain newick> <out.csv> [pretrained_models dir]
import os, sys, re, itertools, pickle, warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"; os.environ["CUDA_VISIBLE_DEVICES"] = ""
warnings.filterwarnings("ignore")
import numpy as np
from ete3 import Tree
from phylodeep.encoding import encode_into_summary_statistics
src, dst = sys.argv[1], sys.argv[2]
rows = []
with open(src) as f:
    next(f)
    for line in f:
        p, nwk = line.rstrip("\n").split("\t")
        k = itertools.count()
        nwk = re.sub(r'([(,]):', lambda m: f'{m.group(1)}t{next(k)}:', nwk)
        tree = Tree(nwk, format=1)
        tree.prune(tree.get_leaves(), preserve_branch_length=True)
        df, resc = encode_into_summary_statistics(tree, float(p))
        rows.append(list(df.iloc[0]) + [resc])
X = np.array(rows, dtype=float)
np.savetxt(dst, X, delimiter=",", fmt="%.17g")
print(X.shape)
if len(sys.argv) > 3:
    d = sys.argv[3]
    sc = pickle.load(open(os.path.join(d, "scalers", "BD_SMALL_FFNN.pkl"), "rb"))
    np.savetxt(dst.replace(".csv", "_scaler.csv"), np.vstack([sc.mean_, sc.scale_]), delimiter=",", fmt="%.17g")
    from tensorflow import keras
    m = keras.models.load_model(os.path.join(d, "models", "BD_SMALL_FFNN.h5"), compile=False)
    out = m.predict(sc.transform(X[:, :-1]), verbose=0)
    np.savetxt(dst.replace(".csv", "_ffnn.csv"), out.astype(float), delimiter=",", fmt="%.10g")
    print("scaler", sc.mean_.shape, "ffnn", out.shape)
