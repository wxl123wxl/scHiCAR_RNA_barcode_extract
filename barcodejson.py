#!/usr/bin/env python3

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--b_dir", help="Required. the FULL path to the fastq folder")
args = parser.parse_args()

assert args.b_dir is not None, "please provide the path to the fastq folder"


## default dictionary is quite useful!

FILES = defaultdict(list)
# collections.defaultdict(list)
## build the dictionary with full path for each fastq.gz file
for root, dirs, files in os.walk(args.b_dir):
    for file in files:
        if file.endswith("_barcode.txt"):
            full_path = join(root, file)
            #R1 will be forward reads, R2 will be reverse reads
            m = re.search(r"(.*)_barcode.txt",file) ## get the sample lane
            if m:
                sample = m.group(1) 
                FILES[sample].append(full_path)


print()
print ("total {}  barcode will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
print ("------------------------------------------")
print("check the barcode.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('barcode.json', 'w').writelines(js)


