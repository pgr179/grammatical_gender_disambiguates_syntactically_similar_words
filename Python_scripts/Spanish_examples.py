#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 11:00:27 2021

@author: pgr
"""


import pandas as pd
import entropies7
import numpy as np

# Load Spanish UDT data
df = pd.read_csv('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-es.csv')

# Assign every word a frequency of 1 so that they sum when we group
df['frequency'] = 1

# Get rid of certain parts of speech
df = df[df.upos == 'NOUN']

# Create gender column if it doesn't already have it (to prevent error)
if "gender" not in df.columns:
    df["gender"] = "none"

# Group by lemmas and parts of speech    
grouped = df.fillna('none').groupby(['lemma','upos'], as_index = False)
lemmadf = grouped.sum()
lemmadf["gender"] = grouped["gender"].agg(lambda x: pd.Series.mode(x)[0])["gender"]

# Check columns
for x in df.columns:
    print(x)

# Verify frequencies match expected
lemmadf.frequency[lemmadf.lemma == 'oro']
lemmadf.frequency[lemmadf.lemma == 'paz']
lemmadf.frequency[lemmadf.lemma == 'medalla']

# Reduce dataframe to just the words we are interested in
words = ['oro', 'paz', 'medalla']
lemmadf = lemmadf.loc[lemmadf['lemma'].isin(words)]
len(lemmadf)

# Get rid of left- and right-specific syntactic relations
lemmadf.drop(list(lemmadf.filter(regex = 'left')), axis = 1, inplace = True)
lemmadf.drop(list(lemmadf.filter(regex = 'right')), axis = 1, inplace = True)
lemmadf.drop(list(lemmadf.filter(regex = 'Unnamed')), axis = 1, inplace = True)

# Cast syntactic information into a list and make that list a new column of the dataframe
synvec = []
for col in lemmadf:
    if ("dep_" in col or "head_" in col):
        synvec.append(col)
lemmadf['syn'] = lemmadf[synvec].values.tolist()

# Apply entropy corrections and turn frequencies into probability distribution
lemmadf.syn = lemmadf.syn.map(lambda x: entropies7.FreqShrink(x))

#Save to file
lemmadf.to_csv('/Users/pgr/Documents/Dissertation/Python_scripts/Span_examples.csv')
