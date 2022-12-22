# -*- coding: utf-8 -*-
"""
Part I:
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

#the parameter csvfile represent the hit table in csv format. If you downloaded it, then
#it does not contain the column names.
# This function set the names of the columns of each blast 
# output in the hit table in csv format according to their respective names.
#If the file has already the respective names, then does nothing.
def inportingHitTable(csvfile):
    header=["query.acc.ver","subject.acc.ver","%.identity","alignment.length","mismatches",
            "gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score","%.positives"]
    with open(csvfile, 'r',newline='') as infile:
        reader = csv.reader(infile)
        if collections.Counter(next(reader)) != collections.Counter(header):
            #temp.* is a temporary file
            with open('temp.csv', 'w',newline='') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(header)
                writer.writerows(reader)
                infile.close()
                outfile.close()
                os.remove(csvfile)
                os.renames("temp.csv", csvfile)
    return 
#This function is a preprocess of BBH methodology
# It returns a join dataframe from the two blast results
#in a comparable way and define as best Hits under the criteria mention below.
#The best hits for each query protein are return for searching best bidirectional hits.
#respective to minimum e-value
def bestJoinBlastHit(file1,file2):
    hits1 = pd.read_csv(file1)
    hits2= pd.read_csv(file2)
    #following line based on:https://stackoverflow.com/questions/61172737/how-to-get-initial-rows-indexes-from-df-groupby
    minEvaluehits1=hits1.groupby(["query.acc.ver"])["evalue"].transform(min)
    BestHits1=hits1.loc[hits1['evalue']==minEvaluehits1]
    minEvaluehits2=hits2.groupby(["query.acc.ver"])["evalue"].transform(min)
    BestHits2=hits2.loc[hits2['evalue']==minEvaluehits2]
    swapBestHits2=swap_columns(BestHits2,'query.acc.ver','subject.acc.ver')
    swapBestHits2.rename(columns={swapBestHits2.columns[0]:"query.acc.ver",swapBestHits2.columns[1]:"subject.acc.ver"},inplace=True)
    bestHits1TwoColumns=BestHits1[['query.acc.ver','subject.acc.ver']]
    swapBestHits2TwoColumns=swapBestHits2[['query.acc.ver','subject.acc.ver']]
    doubleHits=pd.concat([bestHits1TwoColumns,swapBestHits2TwoColumns])
    return doubleHits

    
    
#This function is a preprocess of gene-based analysis methodology
# It returns a join dataframe from the two blast results
#in a comparable way and define best Hits that complies with  the threshold  mention below.
# # Quality score filter parameters:
# Rule of thumb: bit score>= 50-indicator sequence similarity
#"e-value <= 0.01-indicator of number of expected hits of similar quality (score) that could be found just by chance."
#(https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html)
#Under this consideration, all blast results that suceed this threshold are consider
#good hits for sequence homology (including best hits).
#The reason I chose 0.01 as the evalue threshold is that there literature suggesting 
#that evalue smaller to 0.01 can start to be consider to had some homology between them.
#However, weakening this criteria to its limit allows also to have a greater chance to find double hits 
#Reference:
# https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html
def JoinBlastHit(file1,file2):
    hits1 = pd.read_csv(file1)
    hits2= pd.read_csv(file2) 
    bestHits1=(hits1[(hits1["evalue"]<=0.01) & (hits1["bit.score"]>=50)])
    bestHits2=(hits2[(hits2["evalue"]<=0.01) & (hits2["bit.score"]>=50)])
    swapBestHits2=swap_columns(bestHits2,'query.acc.ver','subject.acc.ver')
    swapBestHits2.rename(columns={swapBestHits2.columns[0]:"query.acc.ver",swapBestHits2.columns[1]:"subject.acc.ver"},inplace=True)
    bestHits1TwoColumns=bestHits1[['query.acc.ver','subject.acc.ver']]
    swapBestHits2TwoColumns=swapBestHits2[['query.acc.ver','subject.acc.ver']]
    doubleHits=pd.concat([bestHits1TwoColumns,swapBestHits2TwoColumns])
    return doubleHits
#This method search and filters for the results of the best bidirectional hits inside the the gene-base analysis.
#To further be able to propose a corthologous pair that also is has a bidirectional best hit.
def bestHitsCoorthologyCandidates(bestdoubleHits,doubleHits):
    #this identified all those bbh inside gene base analysis output
    searchBBH = pd.merge(bestdoubleHits, doubleHits, how='left', indicator='exists')
    Candidates1=searchBBH[['query.acc.ver','subject.acc.ver']].loc[searchBBH['exists']=="both"]
    #this retrieve the bbh pair but also the multiple targets to its query( based on gene base analysis output)
    collectionBBH=pd.merge(doubleHits,Candidates1['query.acc.ver'],how='left',indicator='exists')
    Candidates2=collectionBBH[['query.acc.ver','subject.acc.ver']].loc[collectionBBH['exists']=="both"]
    return Candidates2
#Receives the (non)best join hits dataframe of both blast files from 
# the bestJoinBlastHit or JoinBlastHit functions
#This function looks for repeats of query and target if found then exist a bidirectional hit.
def bestBidirectionalHits(doubleHits):
    #bbh:best bidirectional hits
    duplicatesbbh=doubleHits[doubleHits.duplicated(keep='first')]
    bbh=duplicatesbbh.drop_duplicates(keep='first')
    return  bbh
    """
    As this species are not closely related. If we find more that one bidirectional hit from a query protein A
    to multiple target proteins t1,t2,t3,tn such vector called T belonging to species B are likely to be more similar to each other
    than to the query protein from species A. Hence from the definition of slide 17 session 2, the results from this function
    gives you the query protein that are candidate for coorthology. Because targets of this query are likely to be paralogs and the relation
    by BBH method is by orthology.
    """
def coorthologues(Candidates):
    #this returns the rows for which the queries has more than one target
    coortholog=Candidates.loc[Candidates.iloc[:,0].duplicated(keep=False)]
    #All this query proteins are involve in coorthology
    uniqueQuery=coortholog.iloc[:,0].unique()    
    return (coortholog,uniqueQuery)

#This function swaps two columns in a pandas dataframe.
#Used on BestJoinBlastHit function.
#extracted from: https://www.statology.org/swap-columns-pandas/
def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df
"""
In between is obviously part II, which does not need coding.However, part III
follows the continuity of the previous part and starts working over the multiple alignment file 
generated in part II for the species tree. 
"""
#Part III
#Basic method that reads clustal alignment (*.aln) and returns as a 2D array that contains
#the id and the sequence of each species
def reading_alignment(multipleAlignment_aln):
    with open(multipleAlignment_aln, 'r',encoding="utf-8") as file:
        lines=[]
        for line in file:
        #grep tha split by column considering variable spaces
            l1=re.split('\s+',line.strip())
            lines.append(l1)
        #In general the following lines just filter the multiple alignment file in order to
        #keep the ids and the respective sequence of the id
        lines= [ele for ele in lines if len(ele) == 2]
        i=0
        while i <=len(lines)-1:
            #grep pattern based on 
            #https://stackoverflow.com/questions/1507948/search-with-grep-for-words-beginning-and-ending-with
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

#Used in function reading_alignment
#multiple alignment file separates the same sequences many times when 
#comparing with other sequences. This is about to get rid on the read_alignment
#function. However, my code design requires to identify those unique sequences id along
# the file.
#returns an array of unique ids
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

#This is a preprocess step that return a numpy array 
#of the 2D array alignment for the further computations of the sequence conservation by
#Shannon entropy.
#the numpy object does only contain the sequence alignment (the ids are not used)
def numpyAlignment(alignment):
    aligmentArray=[]
    for sequence in alignment:
        aminoAcids=list(sequence[1])
        asciiAmino=[]
        for amino in aminoAcids:
            #ascii coding allows to work with numbers instead of characters
            #facilitates numpy
            asciiAmino.append(ord(amino))
        aligmentArray.append(asciiAmino)
    numpyAlignment=np.array(aligmentArray)
    return numpyAlignment
#This methods is based on Shannon Entropy from information theory
#References:
#https://ieeexplore-ieee-org.kuleuven.e-bronnen.be/document/6773024/authors#authors
#Slides, session 4, slide 7
#Computations of the entropy is based in the dot product from Linear Algebra course.
#transpose are often used because compuational reasons (array storage is easily manipulated on the rows rather than
# the columns).
#returns an array of entropies according to their position
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
        #just a function for applying to a vector
        logTwo=np.vectorize(ma.log2)
        logfrequencies=logTwo(frequencies)
        #entropy dot product
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
        
        
    
    
 #Generate an entropy plot along the sequence for detecting conservation
 #Comparable to the one generated in the webtool in the entropyValidation function.   
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

#This method writes the alignment in a parse format into
#a new file for validating entropy calculations in a  web tool:
#https://www.hiv.lanl.gov/content/sequence/ENTROPY/entropy_one.html
def entropyValidation(alignment):
    with open('EntropyAlignValid.txt','w',encoding="utf-8") as file:
        for sequence in alignment:
            file.write(">"+sequence[0]+'\n')
            file.write(sequence[1]+'\n')    
    return

def main():
    #Part I: Best bidirectional hits
    #One side
    inportingHitTable('acido-staphy.csv')
    
    #Other side
    inportingHitTable('staphy-acido.csv')
    #filtering hits
    geneBaseDoubleHits=JoinBlastHit('staphy-acido.csv','acido-staphy.csv')
    bestDoubleHits=bestJoinBlastHit('staphy-acido.csv','acido-staphy.csv')
    #Best bidirectional hits
    gene_base=bestBidirectionalHits(geneBaseDoubleHits)
    print("These are the best bidirectional hits:")
    bbh=bestBidirectionalHits(bestDoubleHits)
    print(bbh)
    Candidates=bestHitsCoorthologyCandidates(bbh,gene_base)
    # the sample chosen
    
    
    #Coorthologs results
    [coorth,uniqueQuery]=coorthologues(Candidates)
    print("coorthologue results:")
    print(coorth)
    print("list of involve query proteins in coorthology")
    print(uniqueQuery)
    print("this are the sample chosen:")
    print("best bidirectional hit:")
    print(bbh.loc[bbh['query.acc.ver']=="WP_011838606.1"])
    print("coorthologues:")
    print(coorth.loc[coorth['query.acc.ver']=="WP_011838606.1"])
    #Part III: sequence conservation.
    alignment=reading_alignment('multiplealignment.aln')
    npalignment=numpyAlignment(alignment)
    #entropy  calculations
    entropies=entropy(npalignment)
    entropiesArray=dicEntropies(npalignment,entropies)
    print("entropy by aminoacids (first sequence as reference):")
    print(entropiesArray)
    print("generating plot")
    print("printing...")
    entropyPlot(entropies)
    print("your figure has been exported")
    print("generating validation file...")
    entropyValidation(alignment)
    print("validation file has been exported")

    return
if __name__ == "__main__":
    main()
"""
References:
    Best bidirectional hit and gene-based analysis (in Blast section) definition employed 
    http://gencolors.leibniz-fli.de/documentation.html
    Bits core and e value threshold for gene based analysis
    https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html
    https://www.metagenomics.wiki/tools/blast/evalue
    Numpy docs
    https://numpy.org/doc/
    Matplot doc
    https://matplotlib.org/3.6.2/index.html
    Python doc
    https://docs.python.org/3/index.html
    Pandas doc
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
    
"""
