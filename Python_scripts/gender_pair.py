#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 15:47:53 2021

@author: pgr
"""

# Import packages we need
import glob
import os
import pickle
import Levenshtein
import numpy as np
import scipy
import entropies7
import pandas as pd


# Loop through fastText files processing them one by onea
i = 1
for filepath in glob.glob('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_fastText/*.pickle'):
    with open(filepath, "rb") as f:
        vec_df = pickle.load(f)

    # Isolate language code for use in dataframe
    filename = os.path.basename(filepath)
    langcode = filename.split("-")[1].split(".")[0]
    print("Starting: " + str(langcode))
    
    # Subset to just nouns
    vec_df = vec_df[vec_df.upos == 'NOUN']

    # Get corresponding aggregated file
    with open('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_aggregated_lemma/UDT_aggregated_lemma-{}.pickle'.format(langcode), "rb") as f:
        df = pickle.load(f)
    
    # Restrict to length of at least 10
    df = df[df.frequency >= 10]
    
    # Subset to just nouns and get a list of unique genders
    df = df[df.upos == 'NOUN']
    df = df[df.gender != 'none']
    
    genders = list(set(list(df.gender)))
    
    # If there aren't at least two genders, skip this language
    if len(genders) <= 1:
        i += 1
        continue
    
    print("Length of vector dataframe is: " + str(len(vec_df)))
    print("Length of UDT dataframe is: " + str(len(df)))
    
    # Define a new variable that is the length of the word
    #df['length'] = df['lemma'].apply(lambda x: len(x))
    
    # Merge fastText and UDT data
    df = pd.merge(df, vec_df, how = 'inner', on = ['lemma', 'upos'])
    print("Merged length: " + str(len(df)))
    
    # Apply James-Stein shrinkage to get smoothed probability vectors
    df.syn = df.syn.map(lambda x: entropies7.FreqShrink(x))
    df.syn_left = df.syn_left.map(lambda x: entropies7.FreqShrink(x))
    df.syn_right = df.syn_right.map(lambda x: entropies7.FreqShrink(x))
    df.syn_dep = df.syn_dep.map(lambda x: entropies7.FreqShrink(x))
    df.syn_head = df.syn_head.map(lambda x: entropies7.FreqShrink(x))
    
    # Reset index so that I can use indices below
    df = df.reset_index()
    
    # Rename frequency
    df['frequency'] = df.frequency_x
    
    # Get rid of columns we don't need
    df = df[['index', 'lemma', 'gender', 'frequency', 'dispersion', 'syn', 'syn_left', 'syn_right', 'syn_dep', 'syn_head', 'fastText']]
    
    # Get every possible pairing of lemmas
    df['lang'] = langcode
    pairs = pd.merge(df, df, how = 'inner', on = 'lang')
    
    # Get rid of duplicates and self-comparisons
    pairs = pairs[pairs.index_x > pairs.index_y]
    print("Done merging into pairs dataframe...")
    
    # Do calculations on each pair
    pairs['lev_dist'] = list(map(Levenshtein.distance, pairs['lemma_x'], pairs['lemma_y']))
    print("Done with Levenshtein distance...")
    pairs['syn_dist'] = scipy.spatial.distance.jensenshannon(np.array(list(pairs.syn_x)).transpose(), np.array(list(pairs.syn_y)).transpose())
    pairs.drop(['syn_x', 'syn_y'], axis = 1)
    print("Done with first syntactic distance...")
    pairs['syn_left_dist'] = scipy.spatial.distance.jensenshannon(np.array(list(pairs.syn_left_x)).transpose(), np.array(list(pairs.syn_left_y)).transpose())
    pairs.drop(['syn_left_x', 'syn_left_y'], axis = 1)
    print("Done with second syntactic distance...")
    pairs['syn_right_dist'] = scipy.spatial.distance.jensenshannon(np.array(list(pairs.syn_right_x)).transpose(), np.array(list(pairs.syn_right_y)).transpose())
    pairs.drop(['syn_right_x', 'syn_right_y'], axis = 1)
    pairs['syn_dep_dist'] = scipy.spatial.distance.jensenshannon(np.array(list(pairs.syn_dep_x)).transpose(), np.array(list(pairs.syn_dep_y)).transpose())
    pairs.drop(['syn_dep_x', 'syn_dep_y'], axis = 1)
    pairs['syn_head_dist'] = scipy.spatial.distance.jensenshannon(np.array(list(pairs.syn_head_x)).transpose(), np.array(list(pairs.syn_head_y)).transpose())
    pairs.drop(['syn_head_x', 'syn_head_y'], axis = 1)
    print("Done with syntactic distances...")
    pairs['sem_dist'] = list(map(scipy.spatial.distance.cosine, pairs.fastText_x, pairs.fastText_y))
    pairs.drop(['fastText_x', 'fastText_y'], axis = 1)
    print("Done with semantic distance...")
    
    # Get rid of columns we don't need
    pairs = pairs[['lang', 'lemma_x', 'gender_x', 'frequency_x', 'dispersion_x', 'lemma_y', 'gender_y', 'frequency_y', 'dispersion_y', 'lev_dist', 'syn_dist', 'syn_left_dist', 'syn_right_dist', 'syn_dep_dist', 'syn_head_dist', 'sem_dist']]
    
    # Get value of same gender or different gender
    pairs['gender_sameness'] = (pairs['gender_x']==pairs['gender_y']).astype(int)
    
    # Get per million freqs
    ext = pd.read_csv('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv'.format(langcode))
    l = len(ext)
    pairs['p_m_freq_x'] = pairs.frequency_x*(1000000/l)
    pairs['p_m_freq_y'] = pairs.frequency_y*(1000000/l)    
    
    # Save file
    pairs.to_csv('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_gender_paired/UDT_gender_paired-{}.csv'.format(langcode))
    
    # Progress report
    print('Done processing ' + str(langcode) + ' file...')
    print(str(i) + ' out of 61 languages...')
    print('Length was ' + str(len(pairs)))
    i += 1
