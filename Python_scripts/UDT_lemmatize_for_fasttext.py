#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:59:52 2021

@author: pgr
"""

import glob
import os
import json

# Load list of files without reliable lemma information
with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_files_wo_lemmas.txt", "r") as f:
    bad_files = json.loads(f.read())
  
# Load list of shared language codes
with open("/Users/pgr/Documents/Dissertation/Python_scripts/matches.txt", "r") as f:
    matches = json.loads(f.read())

# Loop through all UDT conllu files
for filepath in glob.glob('/Users/pgr/Documents/Corpora_&_Resources/Universal_Dependencies/ud-treebanks-v2.8/**/*.conllu', recursive=True):
    filename = os.path.basename(filepath)
    langcode = filename.split("_")[0]
    
    # Break if the file doesn't have reliable lemma information, skip it
    if filename in bad_files:
        continue
    
    # Break if language codes is not in matches
    if langcode not in matches:
        continue
    
    # Define a place to put lemmas
    all_words = []
    
    # Collect data
    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if(len(line) > 0 and line[0] != '#'):
                tokens = line.split()
                if ("-" not in tokens[0] and tokens[0].isdigit() and tokens[3] != 'PUNCT'):
                    word = {
                        "form": tokens[1],
                        "lemma": tokens[2],
                        "upos": tokens[3],}
                    all_words.append(word)
                    
    # Check if language file already exists to combine with new lemmas
    if os.path.exists("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized_for_fasttext/UDT_lemmatized_for_fasttext-{}.txt".format(langcode)) == True:
        print("Attaching to existing file...")
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized_for_fasttext/UDT_lemmatized_for_fasttext-{}.txt".format(langcode), "r") as f:
            old = json.loads(f.read())
        new = old + all_words
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized_for_fasttext/UDT_lemmatized_for_fasttext-{}.txt".format(langcode), "w") as f:
            f.write(json.dumps(new))    
    
    # Otherwise create a new file
    else:
        print("Creating new file...")
        with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized_for_fasttext/UDT_lemmatized_for_fasttext-{}.txt".format(langcode), "w") as f:
            f.write(json.dumps(all_words))
    print("Done with another file...")
    
