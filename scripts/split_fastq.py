import sys
import csv
import argparse
import textwrap
import pysam


parser = argparse.ArgumentParser( description='Make regions list for splitting jobs.')
parser.add_argument('--input-fastq', type=str, required=True)
parser.add_argument('--splits', type=int, required=True)
args = parser.parse_args()

with pysam.FastxFile(args.input_fastq) as fh:
    for entry in fh:
        # split into n splits
        n = int(len(entry.sequence)/args.splits)
        if n == 0 or len(entry.sequence) == 0:
            continue
        splits = textwrap.wrap(entry.sequence, n)
        i = 0
        #print(splits)
        for seq in splits:
            if i == args.splits:
                continue
            print("@"+entry.name+"-"+str(i))
            print(seq)
            print("+")
            print("="*len(seq))
            i += 1

