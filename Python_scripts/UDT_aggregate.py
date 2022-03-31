#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads in extracted UDT files and aggregates them by lemmas and parts of speech or by wordforms

@author: pgr
"""

import pandas as pd
import glob
import os
import pickle
import json

# Loop through UDT_extracted files
for filepath in glob.glob('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/*.csv'):
    df = pd.read_csv(filepath)
    filename = os.path.basename(filepath)
    langcode = filename.split("-")[1].split('.csv')[0]    
    
    # Assign every word a frequency of 1 so that they sum when we group
    df['frequency'] = 1
    
    # Get rid of certain parts of speech
    poslist = ['PUNCT', 'SYM', 'X']
    df = df[~df.upos.isin(poslist)]
    
    # Create gender column if it doesn't already have it (to prevent error)
    if "gender" not in df.columns:
        df["gender"] = "none"

# Group by lemmas and parts of speech    
    grouped = df.fillna('none').groupby(['lemma','upos'], as_index = False)
    lemmadf = grouped.sum()
    lemmadf["gender"] = grouped["gender"].agg(lambda x: pd.Series.mode(x)[0])["gender"]
    
    # Set a cutoff for frequency and exclude rows that don't meet the cutoff
    lemmadf = lemmadf[lemmadf.frequency >= 5]
    
    # Define new columns representing syntactic relations or subsets
    synvec = []
    for col in lemmadf:
        if ("dep_" in col or "head_" in col) and ("left" not in col) and ("right" not in col):
            synvec.append(col)
    lemmadf['syn'] = lemmadf[synvec].values.tolist()

    synvec = []
    for col in lemmadf:
        if ("dep_" in col or "head_" in col) and ("left" in col):
            synvec.append(col)
    lemmadf['syn_left'] = lemmadf[synvec].values.tolist()

    synvec = []
    for col in lemmadf:
        if ("dep_" in col or "head_" in col) and ("right" in col):
            synvec.append(col)
    lemmadf['syn_right'] = lemmadf[synvec].values.tolist()

    synvec = []
    for col in lemmadf:
        if ("dep_" in col) and ("left" not in col) and ("right" not in col):
            synvec.append(col)
    lemmadf['syn_dep'] = lemmadf[synvec].values.tolist()

    synvec = []
    for col in lemmadf:
        if ("head_" in col) and ("left" not in col) and ("right" not in col):
            synvec.append(col)
    lemmadf['syn_head'] = lemmadf[synvec].values.tolist()
    
    # Remove unwanted columns
    lemmadf = lemmadf[['lemma', 'upos', 'frequency', 'gender', 'syn', 'syn_left', 'syn_right', 'syn_dep', 'syn_head']]
    
    # Load and transform semantic vectors
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_word2vec/UDT_word2vec-{}.pickle".format(langcode), "rb") as f:
        sem_vectors = pickle.load(f)
    vector_df = pd.DataFrame.from_dict(sem_vectors, orient = "index")
    vector_df['sem1'] = vector_df.values.tolist()
    vector_df['lemma'] = vector_df.index.str.rsplit('_', 1).str[0]
    vector_df['upos'] = vector_df.index.str.rsplit('_', 1).str[1]
    vector_df = vector_df[['lemma', 'upos', 'sem1']]
    
    if (len(lemmadf)) == (len(vector_df)):
        print('MATCH')
    else:
        print('CAUTION!!!!')
        print('Length of lemma is ' + str(len(lemmadf)))
        print('Length of vector is ' + str(len(vector_df)))
    
    # Load and transform dispersions
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_dispersions/UDT_dispersions-{}.txt".format(langcode), "r") as f:
        dispersions = json.loads(f.read())
    disp_df = pd.DataFrame.from_dict(dispersions, orient = "index")
    disp_df['dispersion'] = disp_df[disp_df.columns[[0]]]
    disp_df['lemma'] = disp_df.index.str.rsplit('_', 1).str[0]
    disp_df['upos'] = disp_df.index.str.rsplit('_', 1).str[1]
    disp_df = disp_df[['lemma', 'upos', 'dispersion']]
    
    if (len(lemmadf)) == (len(disp_df)):
        print('MATCH')
    else:
        print('CAUTION!!!!')
        print('Length of lemma is ' + str(len(lemmadf)))
        print('Length of disp is ' + str(len(disp_df)))

    # Perform merges
    lemmadf = pd.merge(lemmadf, vector_df, how = 'inner', on = ['lemma', 'upos'])    
    lemmadf = pd.merge(lemmadf, disp_df, how = 'inner', on = ['lemma', 'upos'])
    
    # Save aggregated dataframe
    with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_aggregated_lemma/UDT_aggregated_lemma-{}.pickle".format(langcode), "wb") as f:
        pickle.dump(lemmadf, f)
    
    print("Done with another file: " + str(langcode))
