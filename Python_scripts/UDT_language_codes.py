#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:22:20 2021

@author: pgr
"""

import glob
import os
import json
import pandas as pd

# Load list of files without reliable lemma information
with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_files_wo_lemmas.txt", "r") as f:
    bad_files = json.loads(f.read())

lang_codes = []
lang_names = []

# Loop through all UDT files
for filepath in glob.glob('/Users/pgr/Documents/Corpora_&_Resources/Universal_Dependencies/ud-treebanks-v2.8/**/*.conllu', recursive=True):
    filename = os.path.basename(filepath)

    # Break if the file doesn't have reliable lemma information, skip it
    if filename in bad_files:
        print(filename)
        continue
    
    # Append language code
    lang_codes.append(filename.split("_")[0])
    
    # Append language name
    lang_names.append(filepath.split("/")[-2].split("-")[0].split("UD_")[1])

# Save lists into dataframe
df = pd.DataFrame({'Languages': lang_names, 'Codes': lang_codes})
df = df.drop_duplicates()

# Save dataframe of languages and codes
df.to_csv('/Users/pgr/Documents/Dissertation/Python_scripts/UDT_language_codes.csv')

# fastText codes
codes = 'af als am an ar arz as ast av az azb ba bar bcl be bg bh bn bo bpy br bs bxr ca cbk ce ceb ckb co cs cv cy da de diq dsb dty dv el eml en eo es et eu fa fi fr frr fy ga gd gl gn gom gu gv he hi hif hr hsb ht hu hy ia id ie ilo io is it ja jbo jv ka kk km kn ko krc ku kv kw ky la lb lez li lmo lo lrc lt lv mai mg mhr min mk ml mn mr mrj ms mt mwl my myv mzn nah nap nds ne new nl nn no oc or os pa pam pfl pl pms pnb ps pt qu rm ro ru rue sa sah sc scn sco sd sh si sk sl so sq sr su sv sw ta te tg th tk tl tr tt tyv ug uk ur uz vec vep vi vls vo wa war wuu xal xmf yi yo yue zh'
codes = codes.split(' ')

# Check whether UDT language matches with available fastText vectors
df['fastText'] = 'no'
for code in codes:
    df.fastText[df.Codes == code] = 'yes'

# Create list of codes that are shared by UDT and fastText
matches = list(df.Codes[df.fastText == 'yes'])
matches

with open("/Users/pgr/Documents/Dissertation/Python_scripts/matches.txt", "w") as f:
    f.write(json.dumps(matches))