#!/usr/bin/env python

# REVIEW ALL INDICATIONS HERE: https://workflowy.com/#/081884c4c505 

#Imports 
import sys
sys.path.insert(0, 'PROGRAMS/pyranges/')

import pyranges as pr
from pyfaidx import Fasta
import pandas as pd
import copy

# Definitions
def extend_orfs(p, fasta_path, cds_id = 'ID', stops = ['TAG', 'TGA', 'TAA'], chunk_size = 900): 

    #Clean the Pyranges object
    z = p.df
    w = z[['Chromosome', 'Start', 'End', 'Strand', 'ID', 'Parent', 'Feature']][z.Feature=='CDS']
    w['Strand'] = w['Strand'].astype("string") #So the Pyragnes object is correctly labelled as stranded
    p = pr.PyRanges(w)

    print("\nThis is p:")
    print(p)

    #Load Sequence Data from a Fasta file
    fs = Fasta(fasta_path) #pyfaidx_fasta object, fasta sequences

    #Extend Sequence Upstream
    pup = p.extend({"5":chunk_size}) #pyranges upstream (5') 
    print("\nThis is FULL pup:")
    print(pup)

    #Get the Sequence of the Upstream Extension
    pup = pup.subsequence(0,chunk_size) #Subselect as 'End' what was 'Start' on the full sequence, so we only keep the extended fragment
    print("\nThis is SUBSEQUENCED pup:")
    print(pup)

    #Clip Out of Bounds Intervals
    d_c = {k: len(fs[k]) for k in pup.Chromosome.drop_duplicates()} #Dictionary of Chromosome lengths
    print("This is d_c:\n",d_c)

    pup = pr.genomicfeatures.genome_bounds(pup, d_c, clip=True) #clip=True because we want to keep those elements that have one end within the chromosome
    print("\nThis is OUT OF BOUNDS REMOVED pup:")
    print(pup) 

    # get nucleotide sequence per CDS group in pup (see pyranges  spliced_subsequence), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ASK ABOUT THIS AT SOME POINT
    # let's call it ext5 column; call the resulting dataframe, with one row per CDS group, dg 
    # cds_id = string (name of column, default= None) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HOW CAN IT BE NONE????
#    pup = pup.spliced_subsequence(0, chunk_size, by=cds_id) ########################################### THIS INCREASES THE NUMBER OF ROWS, FOR BOTH subsequence() and spliced_subsequence()
#    print("\nThis is SPLICED SUBSEQUENCE & CDS_ID GROUPED pup:")
#    print(pup)
	
    pup.ext5 = pr.get_sequence(pup, pyfaidx_fasta=fs) #USE THE FUNCTION TRANSCRIPT SEQUENCE (get_transcript_sequence in get_fasta in pyranges)
    print("\nThis should contain the pup sequences:")
    print(pup)

   #Get Index of Rightmost Stop
    pup.up_stop_pos = pup.ext5.apply(lambda x: find_pattern(x, stops, position='last'))  
    print("\nThis should contain the INDEX OF THE RIGHTMOST IN FRAME STOP:")
    print(pup[["Chromosome","up_stop_pos"]])

    pup.up_stop_pos = pd.Series(pup.up_stop_pos, dtype='Int32') #So the column data type is int, not float
    
    print("\nIS THE DATA TYPE ISSUE FIXED?")
    print(pup[["Chromosome","up_stop_pos"]])

    print("\nThis should contain N/A values:") #To take a look at the N/A values
    dataframe_p = pup.df
#    print(dataframe_p[dataframe_p["up_stop_pos"].isin([None])])
    print(dataframe_p[dataframe_p["up_stop_pos"].isin([-1])])

    #using pandas dataframe "apply" on dg, get index of leftmost start in ext5 
    #but only to the right of up_stop_pos  --> this is how much to extend up

################################################################################################################################################## APPLY TO DATAFRAME
######################################################################################################################## MAYBE THERE IS NO STOP

   #cut ext5 at the up_stop_pos. If there is no Stop (N/A) keep the whole thing.
   #Generate subsequence starting at rightmost stop codon (if available) 

    #DEFINE NEW SEQUENCE

    #LOOK FOR THE START PATTERN ON SEQUENCE


def find_pattern(seq, patterns, position):
    #Input is sequence and list of patterns
    #Output is the leftmost / rightmost (depending on position) match for any pattern
    #patterns should have the same length #################################################### assert TO CHECK THIS

    if position == 'last': #i.e. rightmost
        position = -1
    elif position == 'first': #i.e. leftmost
        position = 0

    seq_len = len(patterns[0]) #segment length

    l_m = [i for i in range(0,len(seq),seq_len) if seq[i : i + seq_len] in patterns] #list of matches TRY STARTING ON THE RIGHT IF WE WANT THE LAST

    if len(l_m) == 0: #Otherwise it gets out or range if there is no match
        return -1 #So the resulting column can be numerical. It is a sort of NA that does not modify the column datatype
    else:
#        return l_m[position] + 3 ################################################################################# WHAT ABOUT THIS +3?
        return int(l_m[position])   

#Extend Sequence Downstream

#Execution

# Example: SCRIPTS/extend_orfs_joan.py DATA/mini_Drosophila_melanogaster.BDGP6.32.107.gff3 DATA/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa
# ipython: p = pr.read_gff3('DATA/mini_Drosophila_melanogaster.BDGP6.32.107.gff3')
#          fs = Fasta('DATA/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa')

if __name__ == "__main__": 
    arg = sys.argv
    if len(arg) == 3:
        p = pr.read_gff3(arg[1]) #pyranges object       
        extend_orfs(p,arg[2])
    else:
        print("\nPlease, provide 3 arguments.\n")




