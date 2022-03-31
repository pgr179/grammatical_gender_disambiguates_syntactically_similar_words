#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 12:30:04 2021

@author: pgr
"""

import glob
import os
import json
from collections import Counter

bad_files = []
for filepath in glob.glob('/Users/pgr/Documents/Corpora_&_Resources/Universal_Dependencies/ud-treebanks-v2.8/**/*.conllu', recursive=True):
    lemmas = []
    with open(filepath, "r", encoding="utf-8") as file:
        filename = os.path.basename(filepath)
        langcode = filename.split("_")[0]
        for line in file:
            line = line.strip()
            if(len(line) > 0 and line[0] != '#'):
                tokens = line.split()
                if ("-" not in tokens[0] and tokens[0].isdigit() and tokens[3] != 'PUNCT'):
                    lemma = tokens[2]
                    lemmas.append(lemma)
    c = Counter(lemmas).most_common(5)
    for tup in c:
        if "_" in tup:
            bad_files.append(filename)

with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_files_wo_lemmas.txt", "w") as f:
    f.write(json.dumps(bad_files))
