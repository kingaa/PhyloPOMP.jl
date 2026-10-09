# Reference only (Python is never the pipeline path): phylodeep's CBLV ("most recent") encoding for each tree, as encode_generated_trees.py does it.
# usage: julia --project=pipeline -e "...newick(g; extended=false)..." to make <in.tsv> (columns p_sample, plain newick), then
#        venv/bin/python -I scripts/phylodeep_cblv_oracle.py <in.tsv> <out.csv>   (out: CBLV columns, then the rescale factor)
import sys
import re
import itertools
import numpy as np
from ete3 import Tree
from phylodeep.encoding import encode_into_most_recent

src, dst = sys.argv[1], sys.argv[2]
rows = []
with open(src) as f:
    next(f)
    for line in f:
        p, nwk = line.rstrip("\n").split("\t")
        k = itertools.count()
        nwk = re.sub(r'([(,]):', lambda m: f'{m.group(1)}t{next(k)}:', nwk)   # ete3 needs leaf names
        tree = Tree(nwk, format=1)
        tree.prune(tree.get_leaves(), preserve_branch_length=True)
        df, resc = encode_into_most_recent(tree, float(p))
        rows.append(list(df.iloc[0]) + [resc])
np.savetxt(dst, np.array(rows, dtype=float), delimiter=",", fmt="%.15g")
print(len(rows), "trees,", len(rows[0]), "columns")
