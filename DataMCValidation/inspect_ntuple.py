"""Print the branch layout and value ranges of the converted OPAL ntuples.

Used once to write NtupleVariables.md; kept so the numbers can be reproduced.
    python3 inspect_ntuple.py [events_per_tree]
"""
import sys

import awkward as ak
import numpy as np
import uproot

MERGED = "/raid5/data/yjlee/OPAL/converted/1994/merged/"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 50000


def describe(tree, name):
    print(f"\n===== {name}: {tree.num_entries} entries, {len(tree.keys())} branches")
    arrays = tree.arrays(entry_stop=N)
    for key in tree.keys():
        branch = tree[key]
        values = arrays[key]
        interp = str(branch.interpretation)
        if values.ndim == 1:  # scalar per event
            v = ak.to_numpy(values)
            uniq = np.unique(v)
            extra = f"unique={uniq.tolist()}" if len(uniq) <= 8 else f"min={v.min()} max={v.max()} nunique={len(uniq)}"
            if v.dtype.kind == "f":
                extra += f" mean={v.mean():.4g}"
            print(f"  {key:28s} {interp:34s} {extra}")
        else:  # jagged per particle
            flat = ak.to_numpy(ak.flatten(values))
            if len(flat) == 0:
                print(f"  {key:28s} {interp:34s} (empty)")
                continue
            if flat.dtype.kind in "iub":
                uniq = np.unique(flat)
                extra = f"unique={uniq.tolist()}" if len(uniq) <= 12 else f"min={flat.min()} max={flat.max()} nunique={len(uniq)}"
            else:
                extra = f"min={flat.min():.4g} max={flat.max():.4g} mean={flat.mean():.4g}"
            print(f"  {key:28s} {interp:34s} {extra}")


for filename, trees in [("OPAL_1994_data.root", ["t"]), ("OPAL_1994_mc_jt74mh.root", ["t", "tgen", "tgenBefore"])]:
    f = uproot.open(MERGED + filename)
    print(f"\n########## {filename}: keys = {f.keys()}")
    if "opal_metadata" in f:
        print("opal_metadata:", f["opal_metadata"])
    for t in trees:
        describe(f[t], f"{filename}:{t}")
