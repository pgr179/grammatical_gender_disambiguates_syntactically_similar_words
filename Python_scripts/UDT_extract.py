# -*- coding: utf-8 -*-
"""
This script reads in conllu files from the Universal Dependencies Treebanks and exports a file for each language containing information about each instance of each word, including syntactic information.

Credit: Phillip Rogers & Sherwin Lai
"""

import pandas as pd
import glob
import os
import json

# Load list of files without reliable lemma information
with open("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_files_wo_lemmas.txt", "r") as f:
    bad_files = json.loads(f.read())


# Loop through all UDT files
for filepath in glob.glob('/Users/pgr/Documents/Corpora_&_Resources/Universal_Dependencies/ud-treebanks-v2.8/**/*.conllu', recursive=True):
    filename = os.path.basename(filepath)
    langcode = filename.split("_")[0]

    # Break if the file doesn't have reliable lemma information, skip it
    if filename in bad_files:
        print(filename)
        continue
      
    # Define a list for collected data
    sentences = []
    
    # Collect the data
    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if ("# sent_id" in line):
                sentences.append([])
            elif(len(line) > 0 and line[0] != '#'):
                tokens = line.split()
                if ("-" not in tokens[0] and tokens[0].isdigit()):
                    word = {
                        "lang": langcode,
                        "index": int(tokens[0]),
                        "form": tokens[1],
                        "lemma": tokens[2],
                        "upos": tokens[3],
                        "xpos": tokens[4],
                        "feats": tokens[5].split('|'),
                        "head": -1 if (tokens[6] == '_' or tokens[6] == "NN" or not tokens[6].isdigit()) else int(tokens[6]),
                        "deprel": tokens[7],
                        "deps": tokens[8],
                        "misc": tokens[9].split('|')
                    }
                    for feat in word["feats"]:
                        if "Gender" in feat:
                            word["gender"] = feat.split("=")[1]
                    sentences[len(sentences) - 1].append(word)
    print("Done reading in " + filename + "...")

    allwords = []
    for sentence in sentences:
        for word in sentence:
            if word["head"] > 0:
                key = "dep_" + word["deprel"]
                if 0 < word["head"] <= len(sentence):
                    word2 = sentence[word["head"] - 1]
                    word[key] = 1
                    if word2["index"] - word["index"] > 0:
                        word[key + "_right"] = 1
                    else:
                        word[key + "_left"] = 1
            for otherword in sentence:
                if otherword["head"] == word["index"]:
                    key = "head_" + otherword["deprel"]
                    if key in word:
                        word[key] += 1
                    else:
                        word[key] = 1
                    if otherword["index"] - word["index"] > 0:
                        newkey = key + "_right"
                        if newkey in word:
                            word[newkey] += 1
                        else:
                            word[newkey] = 1
                    else:
                        newkey = key + "_left"
                        if newkey in word:
                            word[newkey] += 1
                        else:
                            word[newkey] = 1
        allwords += sentence
    
    df = pd.DataFrame(allwords)
    for col in df.select_dtypes(include=["float64"]):
        df[col].fillna(0, inplace = True)
        df[col] = df[col].astype("int64")
    
    if os.path.exists("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv".format(langcode)) == True:
        print("Attaching to existing file...")
        df2 = pd.read_csv(("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv".format(langcode)))
        newdf = pd.concat([df, df2])
        for col in newdf.select_dtypes(include=["float64"]):
            newdf[col].fillna(0, inplace = True)
            newdf[col] = newdf[col].astype("int64")
        newdf.to_csv("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv".format(langcode))
    else:
        print("Creating new file...")
        df.to_csv("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_extracted/UDT_extracted-{}.csv".format(langcode))
    
    print("Done reading in " + str(filename) + " of " + str(len(allwords)) + " words")
