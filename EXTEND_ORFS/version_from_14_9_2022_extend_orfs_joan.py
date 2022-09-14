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

########################################################################################################################################################################################################

#### I AM CONFUSED ABOUT THE NEGATIVE STRAND

#    prova_pup = copy.deepcopy(pup)
#    prova_pup.seq = pr.get_sequence(prova_pup, pyfaidx_fasta=fs)

#    dataframe_p = prova_pup.df
#    print("\nI AM TRYING TO UNDERSTAND HOW THE NEGATIVE STRAND OPERATES:")
#    print(dataframe_p[["Chromosome","Start","End","Strand","seq"]])

#    prova_pup = prova_pup.subsequence(0,3)
#    prova_pup.seq = pr.get_sequence(prova_pup, pyfaidx_fasta=fs)

#    dataframe_p = prova_pup.df
#    print("\nSUBSEQUENCED:")
#    print(dataframe_p[["Chromosome","Start","End","Strand","seq"]])


########################################################################################################################################################################################################

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

    print("This should contain N/A values:") #To take a look at the N/A values
    dataframe_p = pup.df
    print(dataframe_p[dataframe_p["up_stop_pos"].isin(["N/A"])])

    #using pandas dataframe "apply" on dg, get index of leftmost start in ext5 
    #but only to the right of up_stop_pos  --> this is how much to extend up

######################################################################################################################## MAYBE THERE IS NO STOP

   #cut ext5 at the up_stop_pos. If there is no Stop (N/A) keep the whole thing.
   #Generate subsequence starting at rightmost stop codon (if available) 

    #Get the Sequence of the Upstream Extension
    #DEAL WITH N/A
#    dataframe_p = pup.df

    print("\nThis should be pup all 32:")
    print(pup)
    pup.old_start = pup.Start
    pup.temporary_start = pup.up_stop_pos.apply(lambda x: 0 if (x == 'N/A') else x) #So it keeps the whole extended sequence if no STOP codon was found ######### PREVIOUS ATTEMPTS BELOW
#    pup.Start = pup.up_stop_pos.apply(lambda x: pup.Start if (x == 'N/A') else x) #This does not work, treats pup.Start as the whole Series
#    pup = pup.subsequence(pup.temporary_start) #SUBSEQUENCE CAN'T USE SERIES :(

    print("\nTHIS SEEMS TO BE STILL ALL 32:")
    print(pup)

    pup.Start = pup.Start + pup.temporary_start ###################################################### A different line for - Strand, pup.End
################################################################################################################################################## APPLY DATAFRAME
    print("\nThis should be pup with some 64:")
    print(pup)

    pup.Start = pd.to_numeric(pup.Start, downcast='integer')
#    pup = pup.subsequence(0) ################################################################################################################################################ DOES THIS MAKE SENSE?

    print("\nIS IT SOLVED???????????????????????????????:")
    print(pup)


    print("\nThis should contain OLD and TEMPORARY start columns (as df):")
    dataframe_p = pup.df
    print(dataframe_p[["Chromosome","Start","End", "up_stop_pos", "old_start","temporary_start"]])   

    print("\nThis should contain N/A and 0s at up_st and temporary_start:") #To take a look at the N/A values
    print(dataframe_p[dataframe_p["temporary_start"].isin([0])][["Chromosome","Start","End", "up_stop_pos", "old_start","temporary_start"]])

    #DEFINE NEW SEQUENCE
    print("\nTake a look at ext5:")
    dataframe_p = pup.df
    print(dataframe_p[["Chromosome","Start","End","Strand","ext5"]])

    print("\nTake a look at ext5 MODIFIED:")
    pup.ext5 = pr.get_sequence(pup, pyfaidx_fasta=fs) #USE THE FUNCTION TRANSCRIPT SEQUENCE (get_transcript_sequence in get_fasta in pyranges)
    dataframe_p = pup.df
    print(dataframe_p[["Chromosome","Start","End","Strand","ext5"]]) ####################################################################################### TROUBLE WITH - STRAND
                                                                     ############################################# I THOUGHT IT COULD BE RELATED TO int64 / int32. BUT IT IS NOT
    print(pup)

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
        return 'N/A' ####################################################################### SUBSTITUTE THIS FOR None CHECK IF THIS IS THEN CONVERTED TO NUMPY N/A
    else:
#        return l_m[position] + 3 ################################################################################# WHAT ABOUT THIS +3?
        return l_m[position]   


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




