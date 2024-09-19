#!/usr/bin/env python3

#By Emre Karakoc, 2020., modified Shriram Bhosle 2024 to use cram file  input
#This is for counting the sgRNA dual CRISPR KO screen reads againts the library
#The reads are aligned to the library with one edit distance and reported the stats
#about the alignment

import time
import sys, os, scipy, numpy
from scipy import stats
from collections import Counter
from Bio.Seq import Seq
import pandas as pd
import pysam

# The Trie data structure keeps a set of words, organized with one node for
# each letter. Each node has a branch for each letter that may follow it in the
# set of words.

library_filename = sys.argv[1]
cram_input_file = sys.argv[2]
myprefix = sys.argv[3]
stat_outputfile = "results/"+myprefix + ".stats"
count_outputfile = "results/"+myprefix + ".counts"
MAX_COST = 0


def intersection(lst1, lst2): 
    temp = set(lst2) 
    lst3 = [value for value in lst1 if value in temp] 
    return lst3 


# Keep some statistics for TRIE structure

statf = open(stat_outputfile, "w")
countf = open(count_outputfile, "w")

parsed_csv = pd.read_csv(library_filename,sep="\t", low_memory=False)

print("The library parameters are read")
print("Start loading the library to the trie")

# read library file into a trie
sgRNA1dic = {}
sgRNA2dic = {}

for index in range(0,len(parsed_csv)):    
    myword = parsed_csv.sgRNA1_WGE_Sequence[index]
    if(len(myword)>19):
        mysgRNA = myword[1:20]
    else:
        mysgRNA = myword[0:19]
    #trie.insert(mysgRNA,parsed_csv.sgRNA1_WGE_ID[index],parsed_csv.ID[index],parsed_csv.sgRNA1_Approved_Symbol[index],"F",1)
    if mysgRNA in sgRNA1dic:
        sgRNA1dic[mysgRNA][1].append(parsed_csv.ID[index])
    else:
        sgRNA1dic[mysgRNA]=[parsed_csv.sgRNA1_WGE_ID[index],[parsed_csv.ID[index]],parsed_csv.sgRNA1_Approved_Symbol[index],"F"]
    seq = Seq(mysgRNA)
    reversecomp = seq.reverse_complement()
    
    #trie.insert(str(reversecomp),parsed_csv.sgRNA1_WGE_ID[index],parsed_csv.ID[index],parsed_csv.sgRNA1_Approved_Symbol[index],"R",1)
    if str(reversecomp) in sgRNA1dic:
        sgRNA1dic[str(reversecomp)][1].append(parsed_csv.ID[index])
    else:
        sgRNA1dic[str(reversecomp)]=[parsed_csv.sgRNA1_WGE_ID[index],[parsed_csv.ID[index]],parsed_csv.sgRNA1_Approved_Symbol[index],"R"]
    myword = parsed_csv.sgRNA2_WGE_Sequence[index]
    if(len(myword)>19):
        mysgRNA = myword[1:20]
    else:
        mysgRNA = myword[0:19]
    #trie.insert(mysgRNA,parsed_csv.sgRNA2_WGE_ID[index],parsed_csv.ID[index],parsed_csv.sgRNA1_Approved_Symbol[index],"R",2)
    if mysgRNA in sgRNA2dic:
        sgRNA2dic[mysgRNA][1].append(parsed_csv.ID[index])
    else:
        sgRNA2dic[mysgRNA]=[parsed_csv.sgRNA2_WGE_ID[index],[parsed_csv.ID[index]],parsed_csv.sgRNA2_Approved_Symbol[index],"F"]
    seq = Seq(mysgRNA)
    reversecomp = seq.reverse_complement()
    #trie.insert(str(reversecomp),parsed_csv.sgRNA2_WGE_ID[index],parsed_csv.ID[index],parsed_csv.sgRNA1_Approved_Symbol[index],"F",2)
    if str(reversecomp) in sgRNA2dic:
        sgRNA2dic[str(reversecomp)][1].append(parsed_csv.ID[index])
    else:
        sgRNA2dic[str(reversecomp)]=[parsed_csv.sgRNA2_WGE_ID[index],[parsed_csv.ID[index]],parsed_csv.sgRNA2_Approved_Symbol[index],"R"]
    
    

# The search function returns a list of all words that are less than the given
# maximum distance from the target word

# This recursive helper is used by the search function above. It assumes that
# the previousRow has been filled in already.


start = time.time()

total_pair_count=0
match=0
matchFF=0
matchRF=0
matchRR=0
swapFF=0
swapFR=0
swapRF=0
swapRR=0
OEM1F=0
OEM1R=0
OEM2F=0
OEM2R=0
nomatch=0

readcount={}

cram = pysam.AlignmentFile(cram_input_file, 'rc', check_sq=False, check_header = False)
reads=cram.fetch(until_eof=True)



# get pairs using iterator next method...

for r1 in reads:
    r2=next(reads)
    #chek if you are using same read pairs...
    if r1.query_name != r2.query_name:
        sys.exit("f{cram_input_file}: Error file is not properly paired")
    # start and get number of based e.g, start from 5 ang get 20 bases
    mytarget1 = r1.query_sequence[0:19]
    mytarget2 = r2.query_sequence[0:19]
    #print(f"{len(mytarget1)}--{len(mytarget2)}")
    total_pair_count=total_pair_count+1
    if mytarget1 in sgRNA1dic:
        if mytarget2 in sgRNA2dic:
            if not set(sgRNA1dic[mytarget1][1]).isdisjoint(sgRNA2dic[mytarget2][1]):
                overlapid = intersection(sgRNA1dic[mytarget1][1],sgRNA2dic[mytarget2][1])[0]
                if sgRNA1dic[mytarget1][3] == "F" and sgRNA2dic[mytarget2][3] == "R":
                    match=match+1
                    if overlapid in readcount:
                        readcount[overlapid]=readcount[overlapid]+1
                    else:
                        readcount[overlapid]=1
                elif sgRNA1dic[mytarget1][3] == "F" and sgRNA2dic[mytarget2][3] == "F":
                    matchFF=matchFF+1
                    new_overlapid = overlapid + "_WO"
                    if new_overlapid in readcount:
                        readcount[new_overlapid]=readcount[overlapid]+1
                    else:
                        readcount[new_overlapid]=1
                elif sgRNA1dic[mytarget1][3] == "R" and sgRNA2dic[mytarget2][3] == "F":
                    new_overlapid = overlapid + "_WO"
                    if new_overlapid in readcount:
                        readcount[new_overlapid]=readcount[overlapid]+1
                    else:
                        readcount[new_overlapid]=1
                    matchRF=matchRF+1
                else:
                    new_overlapid = overlapid + "_WO"
                    if new_overlapid in readcount:
                        readcount[new_overlapid]=readcount[overlapid]+1
                    else:
                        readcount[new_overlapid]=1
                    matchRR=matchRR+1
            else:
                if sgRNA1dic[mytarget1][3] == "F" and sgRNA2dic[mytarget2][3] == "R":
                    swapFR=swapFR+1
                elif sgRNA1dic[mytarget1][3] == "F" and sgRNA2dic[mytarget2][3] == "F":
                    swapFF=swapFF+1
                elif sgRNA1dic[mytarget1][3] == "R" and sgRNA2dic[mytarget2][3] == "F":
                    swapRF=swapRF+1
                else:
                    swapRR=swapRR+1
        else:
            if sgRNA1dic[mytarget1][3] == "F":
                OEM1F=OEM1F+1
            else:
                OEM1R=OEM1R+1
    else:
        if mytarget2 in sgRNA2dic:
            if sgRNA2dic[mytarget2][3] == "F":
                OEM2F=OEM2F+1
            else:
                OEM2R=OEM2R+1
        else:
            nomatch=nomatch+1
        
print("MATCH:",match,total_pair_count,file=statf)
print("MATCHFF:",matchFF,total_pair_count,file=statf)
print("MATCHRF:",matchRF,total_pair_count,file=statf)
print("MATCHRR:",matchRR,total_pair_count,file=statf)
print("SWAPFF:",swapFF,total_pair_count,file=statf)
print("SWAPFR:",swapFR,total_pair_count,file=statf)
print("SWAPRF:",swapRF,total_pair_count,file=statf)
print("SWAPRR:",swapRR,total_pair_count,file=statf)
print("OEM1F:",OEM1F,total_pair_count,file=statf)
print("OEM1R:",OEM1R,total_pair_count,file=statf)
print("OEM2F:",OEM2F,total_pair_count,file=statf)
print("OEM2R:",OEM2R,total_pair_count,file=statf)
print("NOMATCH:",nomatch,total_pair_count,file=statf)

for key, value in readcount.items():
    print(key, value, file=countf)

end = time.time()
print("Search took %g s" % (end - start),file=statf)

statf.close()
countf.close()


