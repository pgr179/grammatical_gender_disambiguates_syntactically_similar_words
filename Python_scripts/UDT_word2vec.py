#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:10:34 2021

@author: pgr
"""

import glob
import os
import json
import pickle
from gensim.models import Word2Vec

# Loop through lemmatized files, saving language code for later
for filepath in glob.glob('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_lemmatized/*.txt'):
    filename = os.path.basename(filepath)
    langcode = filename.split("-")[1].split(".")[0]
    print(langcode)
    
    # Load corpus
    with open(filepath, "r") as f:
        corpus = json.loads(f.read())
        
    # Turn nested lists into straight list of lemmas
    all_lemmas = []
    for i in corpus:
        for j in i:
            all_lemmas.append(j)
    unique_lemmas = list(set(all_lemmas))

    # Train model
    model = Word2Vec(sentences = corpus, vector_size = 100, window = 5, min_count = 1, workers = 7, epochs = 25)
                         
    # Get vectors for each word
    sem_vectors = {}
    for lemma in unique_lemmas:
        sem_vectors[lemma] = model.wv[lemma]
        
    # Save vectors
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_word2vec/UDT_word2vec-{}.pickle".format(langcode), "wb") as f:
        pickle.dump(sem_vectors, f)
    print("Done with another file...")
