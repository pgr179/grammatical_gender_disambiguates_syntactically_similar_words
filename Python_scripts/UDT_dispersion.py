#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:12:41 2021

@author: pgr
"""

import glob
import os
import json
from collections import Counter

# Loop through UDT lemmatized files
for filepath in glob.glob('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/*.txt'):
    with open(filepath, "r") as f:
        corpus = json.loads(f.read())
        filename = os.path.basename(filepath)
        langcode = filename.split("-")[1].split(".")[0]
        print(langcode)

    # Turn nested lists into straight list of lemmas
    lemmalist = []
    for i in corpus:
        for j in i:
            lemmalist.append(j)
    
    # Define length of corpus
    l = len(lemmalist)
    
    # Get counts for each lemma in corpus
    c = Counter(lemmalist)
    
    # Create dictionary of indices
    all_indices = {}
    for index, lemma in enumerate(lemmalist):
        if lemma in all_indices:
            all_indices[lemma].append(index)
        else:
            all_indices[lemma] = [index]
    
    # Get lemmas that occur at least five times
    freq_lemmas = []
    for lemma in c:
        if c[lemma] > 4:
            freq_lemmas.append(lemma)
    len(freq_lemmas)
    
    # Define dictionary for lemmas and dispersions
    dispdict = {}
    
    # Loop through lemmas to calculate ARF dispersion
    for lemma, indices in all_indices.items():
        
        #  frequency of lemma 
        f = c[lemma]
        
        # Get distances (d) between each index, including wrapped distance from last to first
        d = []
        d.append(indices[0] + (l - indices[-1]))
        for index_of_index, index in enumerate(indices):
            if index_of_index > 0:
                d.append(index - indices[index_of_index-1])
        
        # Use ARF dispersion formula
        fp = l/f
        mins = []
        for index, distance in enumerate(d):
            mins.append(min(fp, distance))
        ARF = (f/l) * (sum(mins))
        
        # Add lemma and disperion to dictionary
        dispdict[lemma] = ARF
    
    # Save dispersions
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_dispersions/UDT_dispersions-{}.txt".format(langcode), "w") as f:
        f.write(json.dumps(dispdict))
    print("Done with another file...")