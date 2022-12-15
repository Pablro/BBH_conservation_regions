# -*- coding: utf-8 -*-
"""
Part 1:
Bidirectional Best Hit
"""
import csv 
import os
import collections
import pandas as pd
import re
import numpy as np
import math as ma
import matplotlib.pyplot as plt

#This function set the names of the columns accordingly
# of each bast csv file output.
#if the file has already the respective names, then does anything.
def inportingHitTable(csvfile):
    header=["query.acc.ver","subject.acc.ver","%.identity","alignment.length","mismatches",
            "gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score","%.positives"]
    with open(csvfile, 'r',newline='') as infile:
        reader = csv.reader(infile)
        if collections.Counter(next(reader)) != collections.Counter(header):
            #temp.* is a temporary file
            with open('temp.csv', 'w',newline='') as outfile:
                writer = csv.writer(outfile)
                # create list with old header and new
                writer.writerow(header)
                writer.writerows(reader)
                infile.close()
                outfile.close()
                os.remove(csvfile)
                os.renames("temp.csv", csvfile)
    return 

#This function reads directly csv files 
def readingHitTable(csvfile):
    with open(csvfile, newline='') as f:
        csvreader=csv.reader(f)
        for row in csvreader:
            print(row)
        f.close()
    return

#This function returns a join dataframe from the two blast results
#in a comparable way.
# # Quality score filter parameters:
# Rule of thumb: bit score 50-indicator sequence similiraty
#e value of 0.01-indicator of number of hits by alignment with database
# https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html
def bestJoinBlastHit(file1,file2):
    hits1 = pd.read_csv(file1)
    hits2= pd.read_csv(file2) 
    bestHits1=(hits1[(hits1["evalue"]<=0.01) & (hits1["bit.score"]>=50)])
    bestHits2=(hits2[(hits2["evalue"]<=0.01) & (hits2["bit.score"]>=50)])
    swapBestHits2=swap_columns(bestHits2,'query.acc.ver','subject.acc.ver')
    swapBestHits2.rename(columns={swapBestHits2.columns[0]:"query.acc.ver",swapBestHits2.columns[1]:"subject.acc.ver"},inplace=True)
    #difference in evalues and starting points makes the following methods to not work. However,
    #both balst highlight same inference that both sequences are or not are bidirectional hits and
    #significat according to score filter for homology 
    bestHits1TwoColumns=bestHits1[['query.acc.ver','subject.acc.ver']]
    swapBestHits2TwoColumns=swapBestHits2[['query.acc.ver','subject.acc.ver']]
    doubleHits=pd.concat([bestHits1TwoColumns,swapBestHits2TwoColumns])
    return doubleHits
#Receives the best join hits dataframe of both blast files from 
# the bestJoinBlastHit function
def bestBidirectionalHits(doubleHits):
    
    #bbh=best bidirectional hits
    duplicatesbbh=doubleHits[doubleHits.duplicated(keep='first')]
    bbh=duplicatesbbh.drop_duplicates(keep='first')
    return  bbh
#unique best hits in each blast
def bestNoBidirectionalHits(doubleHits):
    uniqueBestHits=doubleHits.drop_duplicates(keep=False)
    return uniqueBestHits
    """
    As this sepcies are not closely related. If we find more that one bidirectional hit from a query protein A
    to multiple target proteins b1,b2,b3,bn such vector called B belonging to species B are likely to be more similar to each other
    than to the query protein from species A. Hence from the definition of slide 17 session 2, the results from this function
    gives you the query protein that are candidate for coorthology. Because subject query are likely to be paralogs and the relation
    with the BBH is by orthology.
    """
def coorthologues(BBH):
    coortholog=BBH.loc[BBH.iloc[:,0].duplicated(keep=False)]
    #All this query proteins are involve in corthology
    uniqueQuery=coortholog.iloc[:,0].unique()
    return (coortholog,uniqueQuery)
    return
#Summary of findings
def summaryHomologs(bbh,uniqueBestHits):
    countBBH=bbh.shape[0]
    countuniqueBestHits=uniqueBestHits.shape[0]
    counts=[countBBH,countuniqueBestHits]
    report=pd.DataFrame(counts,index=["BBH","Unique"],columns=["Orthologues"])
    print("this is a contingency table of the findings:")
    print(report)
    print("these are Best bidirectional Hits:")
    print(bbh)
    print("these are unique best hits")
    print(uniqueBestHits)
    return
#This function swaps two columns in a pandas dataframe
#extracted from: https://www.statology.org/swap-columns-pandas/
def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df

#Part III
#Simple method that read clustal alignment and returns a 2D array
def reading_alignment(multipleAlignment_aln):
    with open(multipleAlignment_aln, 'r',encoding="utf-8") as file:
        lines=[]
        for line in file:
            l1=re.split('\s+',line.strip())
            lines.append(l1)
        lines= [ele for ele in lines if len(ele) == 2]
        i=0
        while i <=len(lines)-1:
            if re.search("^[a-zA-Z]", lines[i][0]) is None:
                del lines[i]
                i=-1 
            i=i+1
        ids=unique_ids(lines)
        alignment=[]
        for id in ids:
            sequence=""
            for line in lines:
                if id == line[0]:
                    sequence=sequence+line[1].replace(" ","")
            alignment.append([id,sequence])   
    return alignment
#This is based in my code design
#look at function reading_alignment
def unique_ids(alignmentArray):
    idArray=[]
    for sequences in alignmentArray:
        id=sequences[0]
        if len(idArray)!=0:
            if id not in idArray:
                idArray.append(id)
        else:
            idArray.append(id)
    return idArray

#This is a preprocess step that return a numpy array for computationa
#calculations of entropy.
def numpyAlignment(alignment):
    aligmentArray=[]
    for sequence in alignment:
        aminoAcids=list(sequence[1])
        asciiAmino=[]
        for amino in aminoAcids:
            asciiAmino.append(ord(amino))
        aligmentArray.append(asciiAmino)
    numpyAlignment=np.array(aligmentArray)
    return numpyAlignment
#This methods is based on Shannon Entropy from information theory
#https://ieeexplore-ieee-org.kuleuven.e-bronnen.be/document/6773024/authors#authors
#Slides
def entropy(numpyAlignment):
    entropies=[]
    alignmentMatrix=np.matrix(numpyAlignment)
    transposeMtarix=alignmentMatrix.getT()
    numberColumns=transposeMtarix.shape[0]
    for col in range(0,numberColumns):
        column=np.asarray(transposeMtarix[col,:]).flatten()
        columnComp=np.matrix(np.column_stack(np.unique(column,return_counts=True)))
        columnComptrans=columnComp.getT()
        frequencies=columnComptrans[1,:]*(1/np.sum(columnComptrans[1,:]))
        logTwo=np.vectorize(ma.log2)
        logfrequencies=logTwo(frequencies)
        entropy=(-1)*np.dot(frequencies,logfrequencies.getT())
        entropyN=np.asarray(entropy).flatten()[0]
        entropies.append(entropyN)
    return entropies
#This are all the entropies by aminoacids relative to the first sequence.
def dicEntropies(numpyAlignment,entropies):
    entropiesArray=[]
    firstSequence=numpyAlignment[0]
    for amino in range(0,len(firstSequence)):
        entropiesArray.append([chr(firstSequence[amino]),entropies[amino]])
    return entropiesArray
        
        
    
    
 #Entropy plot along the sequence for detecting conservation   
def entropyPlot(entropies):
    
    x=list(range(1,len(entropies)+1))
    y=entropies
    fig, ax = plt.subplots()
    ax.bar(x, y, width=1,color="red", edgecolor="red", linewidth=0.7)
    ax.set(xlim=(0, len(x)), xticks=np.arange(0, len(x),len(x)/6),
           ylim=(0, max(entropies)), yticks=np.arange(0,max(y),max(y)/4 ))
    ax.tick_params(axis="x",labelsize=5)
    ax.tick_params(axis="y",labelsize=5)
    ax.set_ylabel("entropy")
    ax.set_xlabel("sequence")
    ax.set_title("Sequence Conservation")
    fig.tight_layout()
    plt.savefig("entropy-plot.png", format="png")

    
    return 
#This method write the alignment in a parse format into
#a new file for validating entropy calculations in web tool:
#https://www.hiv.lanl.gov/content/sequence/ENTROPY/entropy_one.html
def entropyValidation(alignment):
    with open('EntropyAlignValid.txt','w',encoding="utf-8") as file:
        for sequence in alignment:
            file.write(">"+sequence[0]+'\n')
            file.write(sequence[1]+'\n')    
    return
def main():
    #One side
    inportingHitTable('acido-staphy.csv')
    
    #Other side
    inportingHitTable('staphy-acido.csv')
    #testing hits
    doubleHits=bestJoinBlastHit('staphy-acido.csv','acido-staphy.csv')
    #Best bidirectional hits
    bbh=bestBidirectionalHits(doubleHits)
    
    #Best not bidirectional hits
    uniqueBestHits=bestNoBidirectionalHits(doubleHits)
    #Summary of orthologues hits
    summaryHomologs(bbh,uniqueBestHits)
    #Coorthologs results
    [coorth,uniqueQuery]=coorthologues(bbh)
    print("list of coorthologues")
    print(coorth)
    print("list of involve query proteins in coorthology")
    print(uniqueQuery)
    alignment=reading_alignment('multiplealignment.aln')
    npalignment=numpyAlignment(alignment)
    entropies=entropy(npalignment)
    entropiesArray=dicEntropies(npalignment,entropies)
    print("entropy by aminoacids (first sequence as reference):")
    print(entropiesArray)
    print("generating plot")
    print("printing...")
    entropyPlot(entropies)
    print("your figure has been exported")
    entropyValidation(alignment)

    return
if __name__ == "__main__":
    main()
"""
References:
    http://gencolors.leibniz-fli.de/documentation.html
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101850
    https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html
    https://www.metagenomics.wiki/tools/blast/evalue
    
"""
