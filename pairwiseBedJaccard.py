#! /usr/bin/env python3

# Vivek Rai
# Parker Lab
# vivekrai@umich.edu
#
# Jun 11, 2019

import os
# import logging
import sys
import glob

import docopt
import pybedtools as bt
import pandas

doc = """
    Usage:
        pairwiseBedJaccard (-a <A>)... (-b <B>)... [--options=OPT]

    Options:
    -a <A>          BED files for one set
    -b <B>          BED files for second set
    --options=OPT   Metric to compute (not implemented) [default:jaccard]
"""

opts = docopt.docopt(doc)

# logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")

# logging.info(
#     f"Will compute pairwise Jaccard metrics between {len(opts['-a'])} "
#     f"file(s) in first set and {len(opts['-b'])} in second set."
# )


def scan_expand_args(arg):
    result = []
    for i in arg:
        if "*" in i:
            result.extend(glob.glob(i))
        else:
            result.append(i)
    return result


opts["-a"] = scan_expand_args(opts["-a"])
opts["-b"] = scan_expand_args(opts["-b"])

result = []
for f_a in opts["-a"]:
    f_a_bt = bt.BedTool(f_a)
    for f_b in opts["-b"]:
        jc = f_a_bt.jaccard(bt.BedTool(f_b))
        jc.update({"source": os.path.basename(f_a), "with": os.path.basename(f_b)})
        result.append(jc)


pandas.DataFrame(result).to_csv(sys.stdout, sep="\t", index=False)
