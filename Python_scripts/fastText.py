#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 21:52:23 2021

@author: pgr
"""

import glob
import os
import pandas as pd
import numpy as np
import sys
import tqdm
import pickle


### Load fastText file

for filepath in glob.glob('/Volumes/PGR/FastText_word_vectors/*.txt'):
    filename = os.path.basename(filepath)
    langcode = filename.split(".")[1]

    # Get file size (for later progress bar)
    fileSize = os.path.getsize(filepath)
    
    # No need to preallocate dictionary size
    d = dict()
    
    print("Processing {} ({:0.2f} GB)".format(filepath, fileSize / (10**8)))
    
    with open(filepath, 'r') as file:
        # Create progress bar (negligable impact on processing)
        tq = tqdm.tqdm(total=fileSize, unit_scale=True)
    
        # Skip first line of file & update progress bar
        tq.update(len(file.readline()))
    
        try:
            
            for line in file:
            
                expl = line.split(' ')
        
                # Cast list of strings into np array. Saving a list of strings into the dict
                # consumes all 32GB on my computer, but treating them as an array of floats
                # makes this about 40 MB and takes about 29 sec total on my computer
                d[expl[0]] = np.asarray(expl[1:], dtype=np.float64, order='C')
        
                # Update progress bar
                tq.update(len(line))
        except:
            pass
    # Tidy up progess bar
    tq.close()
    
    print("done; dictionary size {:0.2f} MB".format(sys.getsizeof(d) / (10**9)))
    
    
    ### Load UDT file
    
    # Find appropriate UDT file
    df = pd.read_csv('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv'.format(langcode))   
    
    # Keep only columns that are needed
    df = df[['form', 'lemma', 'upos']]
    
    # Assign every word a frequency of 1 so that they sum when we group
    df['frequency'] = 1
    
    # Aggregate on wordform, lemma, and upos
    df = df.groupby(['form', 'lemma', 'upos'], as_index = False).sum()
    
    # Loop through UDT data to get vectors from hash
    vectors = []
    for word in df.form:
        try:
            vectors.append(d[word])
        except:
            vectors.append(np.NaN)
    df['fastText'] = vectors
    df = df.dropna().reset_index()
    print("Done matching vectors to dataframe...")
    
    # Cast vectors into their own columns
    newdf = pd.DataFrame(df.fastText.tolist())
    
    # Join new columns to original dataframe
    df = pd.concat([df, newdf], axis = 1)
    
    # Multiply vectors by frequency
    df.iloc[:, 6:306] = df.iloc[:, 6:306].mul(df['frequency'], axis = 0)
    
    # Aggregate on lemma and upos
    df = df.groupby(['lemma', 'upos'], as_index = False).sum()

    # Divide vectors by frequency
    df.iloc[:, 4:304] = df.iloc[:, 4:304].div(df['frequency'], axis = 0)
    
    # Coerce vectors back to single array
    df['fastText'] = df.iloc[:, 4:304].values.tolist()
    
    # Get rid of unwanted columns
    df = df[['lemma', 'upos', 'frequency', 'fastText']]
    
    # Restrict to length of at least 10
    df = df[df.frequency >= 10]
    
    # Save file
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_fastText/UDT_fastText-{}.pickle".format(langcode), "wb") as f:
        pickle.dump(df, f)
    print("Done with language: " + str(langcode))
    print("Saved dataframe was length " + str(len(df)))
