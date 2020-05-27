#! /usr/bin/env python3
#
# Vivek Rai
# mail@raivivek.in
# (c) Parker Lab
#
# Jul 8, 2019
# GPLv3


import docopt
import sys
import signal
import pathlib

doc = """
    formatHiCtoCicero

    Prints to STDOUT.

    Usage: formatHiCtoCicero <file>

    Options:
    <file>        input
    -h, --help    Show help
"""


def convert_to_bedpe(file):
    # conns = dict()
    with open(file, "r") as f:
        for line in f:
            loop1, loop2, _ = line.split("\t")
            tab = "\t"
            sys.stdout.write(
                f"{tab.join(loop1.split('_'))}" "\t" f"{tab.join(loop2.split('_'))}\n"
            )


if __name__ == "__main__":
    # Register SIGNIT handler for graceful exit
    # See: https://stackoverflow.com/questions/4205317
    signal.signal(signal.SIGINT, lambda s, f: sys.exit(0))

    opts = docopt.docopt(doc, version="0.1.0")
    file_path = pathlib.Path(opts["<file>"])

    # Exit early if the path does not exist
    if not file_path.exists():
        print("Error: File not found. Exit", file=sys.stderr)
        sys.exit(1)

    convert_to_bedpe(file_path)
