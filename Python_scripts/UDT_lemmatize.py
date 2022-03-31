#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script goes through the UDT treebanks and lemmatizes the corpus by part of speech, creating a list(sentences) of lists (lemmas)

@author: pgr
"""

import glob
import os
import json

# Load list of files without reliable lemma information
with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_files_wo_lemmas.txt", "r") as f:
    bad_files = json.loads(f.read())

# Loop through all UDT conllu files
for filepath in glob.glob('/Users/pgr/Documents/Corpora_&_Resources/Universal_Dependencies/ud-treebanks-v2.8/**/*.conllu', recursive=True):
    filename = os.path.basename(filepath)
    langcode = filename.split("_")[0]
    
    # Break if the file doesn't have reliable lemma information, skip it
    if filename in bad_files:
        print(filename)
        continue
    
    # Define a place to put lemmas
    sentences = []
    
    # Collect data
    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if ("# sent_id" in line):
                sentences.append([])
            elif(len(line) > 0 and line[0] != '#'):
                tokens = line.split()
                if ("-" not in tokens[0] and tokens[0].isdigit() and tokens[3] != 'PUNCT'):
                    lemma = str(tokens[2] + '_' + tokens[3])
                    sentences[len(sentences) - 1].append(lemma)    
    
    # Check if language file already exists to combine with new lemmas
    if os.path.exists("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/UDT_lemmatized-{}.txt".format(langcode)) == True:
        print("Attaching to existing file...")
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/UDT_lemmatized-{}.txt".format(langcode), "r") as f:
            old = json.loads(f.read())
        new = old + sentences
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/UDT_lemmatized-{}.txt".format(langcode), "w") as f:
            f.write(json.dumps(new))    
    
    # Otherwise create a new file
    else:
        print("Creating new file...")
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/UDT_lemmatized-{}.txt".format(langcode), "w") as f:
            f.write(json.dumps(sentences))
    print("Done with another file...")