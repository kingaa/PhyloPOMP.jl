#!/usr/bin/env python3
"""
python_cblv_crosscheck.py

Read SEIR Newick trees from output/crosscheck/trees/,
apply phylodeep CBLV encoding, and write raw (x, y) vectors
(before and after rescaling) to CSV so Julia can compare.

Usage:
    /home/dislam/dl_phylo_env/bin/python scripts/python_cblv_crosscheck.py
"""

import os
import csv
import sys
import numpy as np
from ete3 import Tree

# We replicate the core of phylodeep's encode_into_most_recent but
# capture the raw interleaved vector *before* reshaping/padding,
# and also apply the same rescaling so Julia can compare both.

sys.path.insert(0, os.path.dirname(__file__))

from phylodeep.encoding import rescale_tree, TARGET_AVG_BL

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(SCRIPT_DIR, "..")
CROSSCHECK = os.path.join(BASE_DIR, "output", "crosscheck")
MANIFEST = os.path.join(CROSSCHECK, "manifest.csv")
OUT_DIR = os.path.join(CROSSCHECK, "python_cblv")
os.makedirs(OUT_DIR, exist_ok=True)


def name_tree(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        if not node.is_leaf():
            node.name = f"int_{i}"


def add_dist_to_root(tree):
    tree.dist_to_root = 0
    for node in tree.traverse("preorder"):
        if node.up:
            node.dist_to_root = node.up.dist_to_root + node.dist


def real_polytomies(tre):
    for nod in tre.traverse("postorder"):
        if not nod.is_leaf() and not nod.is_root():
            if nod.dist == 0:
                for child in nod.children:
                    nod.up.add_child(child)
                nod.up.remove_child(nod)


def get_not_visited_anc(leaf):
    while getattr(leaf, "visited", 0) >= len(leaf.children) - 1:
        leaf = leaf.up
        if leaf is None:
            break
    return leaf


def get_deepest_not_visited_tip(anc):
    max_dist = -1
    tip = None
    for leaf in anc:
        if leaf.visited == 0:
            distance_leaf = leaf.dist_to_root - anc.dist_to_root
            if distance_leaf > max_dist:
                max_dist = distance_leaf
                tip = leaf
    return tip


def encode_walk(anc):
    """Generator: interleaved [tip_depth, internal_depth, tip, internal, ...]"""
    leaf = get_deepest_not_visited_tip(anc)
    yield leaf.dist_to_root - anc.dist_to_root  # tip branch
    leaf.visited += 1
    anc = get_not_visited_anc(leaf)
    if anc is None:
        return
    anc.visited += 1
    yield anc.dist_to_root  # internal depth
    yield from encode_walk(anc)


def cblv_raw(nwk_str):
    """
    Returns (x_raw, y_raw, rescale_factor) using the same algorithm
    as phylodeep's encode_into_most_recent, but split into separate
    x (tip) and y (internal) vectors without zero-padding.

    The values are rescaled (mean branch length → 1).
    """
    tree = Tree(nwk_str, format=1)

    # strip single-child root edge (same as phylodeep)
    if len(tree.children) < 2:
        tree = tree.children[0]
        tree.detach()

    real_polytomies(tree)
    rescale_factor = rescale_tree(tree, target_avg_length=TARGET_AVG_BL)

    for node in tree.traverse():
        node.visited = 0
    name_tree(tree)
    add_dist_to_root(tree)

    raw = list(encode_walk(tree))
    # raw is interleaved: [tip0, int0, tip1, int1, ..., tip_last]
    # tip entries at even indices, internal at odd
    x = [raw[i] for i in range(0, len(raw), 2)]
    y = [raw[i] for i in range(1, len(raw), 2)]

    return x, y, rescale_factor


def cblv_unscaled(nwk_str):
    """
    Same walk but WITHOUT rescaling branch lengths.
    Returns (x_raw, y_raw) in the original time units.
    """
    tree = Tree(nwk_str, format=1)

    if len(tree.children) < 2:
        tree = tree.children[0]
        tree.detach()

    real_polytomies(tree)
    # NO rescaling

    for node in tree.traverse():
        node.visited = 0
    name_tree(tree)
    add_dist_to_root(tree)

    raw = list(encode_walk(tree))
    x = [raw[i] for i in range(0, len(raw), 2)]
    y = [raw[i] for i in range(1, len(raw), 2)]

    return x, y


def main():
    with open(MANIFEST) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"Processing {len(rows)} trees from {MANIFEST}")

    summary = []
    for r in rows:
        idx = int(r["index"])
        nwk_path = os.path.join(CROSSCHECK, r["nwk_file"])
        with open(nwk_path) as f:
            nwk = f.read().strip()

        # Unscaled (original time units) — for direct comparison with Julia
        x_unsc, y_unsc = cblv_unscaled(nwk)

        # Rescaled (phylodeep-style)
        x_resc, y_resc, rescale = cblv_raw(nwk)

        # Write unscaled CSV
        csv_path = os.path.join(OUT_DIR, f"tree_{idx:04d}_unscaled.csv")
        with open(csv_path, "w", newline="") as cf:
            w = csv.writer(cf)
            w.writerow(["position", "tip_x", "internal_y"])
            for k in range(len(x_unsc)):
                y_val = y_unsc[k] if k < len(y_unsc) else ""
                w.writerow([k, f"{x_unsc[k]:.12f}", f"{y_val:.12f}" if y_val != "" else ""])

        # Write rescaled CSV
        csv_path_r = os.path.join(OUT_DIR, f"tree_{idx:04d}_rescaled.csv")
        with open(csv_path_r, "w", newline="") as cf:
            w = csv.writer(cf)
            w.writerow(["position", "tip_x", "internal_y"])
            for k in range(len(x_resc)):
                y_val = y_resc[k] if k < len(y_resc) else ""
                w.writerow([k, f"{x_resc[k]:.12f}", f"{y_val:.12f}" if y_val != "" else ""])

        n_tips = len(x_unsc)
        print(f"  tree {idx}: tips={n_tips}  rescale={rescale:.6f}  "
              f"height_unsc={x_unsc[0]:.4f}  height_resc={x_resc[0]:.4f}")

        summary.append({
            "index": idx,
            "n_tips": n_tips,
            "n_internal": len(y_unsc),
            "rescale": rescale,
            "height_unscaled": x_unsc[0],
            "height_rescaled": x_resc[0],
        })

    # Summary CSV
    with open(os.path.join(OUT_DIR, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary[0].keys())
        w.writeheader()
        w.writerows(summary)

    print(f"\nWrote {len(rows)} trees to {OUT_DIR}")


if __name__ == "__main__":
    main()
