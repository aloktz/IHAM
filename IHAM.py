#!/usr/bin/env python
# coding: utf-8

# Author: Alon Katz
# 
# IHAM - Inference of Host Adaptive Mutations
# 
# 
#     

# In[1]:

#imports
from os.path import exists
from Bio import AlignIO, SeqIO
import math
from Bio.Align.AlignInfo import SummaryInfo
import sys
import argparse
import os
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Phylo.Applications import _Fasttree
import subprocess
import json
from ete3 import Tree
import pandas as pd
import time


 #set command line args
parser = argparse.ArgumentParser()
parser.add_argument("-p","--path",help="path to folder of alignment files")
parser.add_argument("-f","--format",help="type the format of your alignment in lower case")
parser.add_argument("-e","--extension",help="type the file extension of your alignment files (case sensitive)")
parser.add_argument("-s","--similarity",help="float value between 0 and 1 - impacts number of sequences sent to FUBAR. Low similarity threshold = more sequences, up to 499 sequences",type=float)

#parse args
args = parser.parse_args()

path=args.path
file_format=args.format
file_extension = args.extension
similarity_threshold = args.similarity
fasta_files=[path+"/"+f for f in os.listdir(path) if f.endswith(file_extension)]

# In[2]:


#getting user input
#file_format = input('what is the format of your alignment file? Type in lower case e.g. "fasta"')
#file_extension = input('what is the extension used by your file format?')
#bool = False
#while bool == False:
#    path = input('please enter the absolute path to your alignment folder')
#    if os.path.isdir(path): 
#        fasta_files=[path+"/"+f for f in os.listdir(path) if f.endswith(file_extension)]
#        bool = True
#        print('file located')
#    else:
#        print('file not found. please try again.')
#
##read_or_parse = input('does your file contain one alignment? type "yes" or "no"')
#print(fasta_files)
#
#for use in read_alignment() and then busted()



# In[3]:


#read alignment
def read_alignment(i: int):
    path = file_data[i][0]+'-'+file_data[i][1]+'-'+file_data[i][2]
    alignment = AlignIO.read(open(path),file_format)
    
    return alignment

#alignment = read_alignment(0)




# In[4]:


#find terminal gaps
def trim_end_gaps(alignment):
    pos_gap = math.inf
    for i in range(len(alignment)):
        s=alignment[i].seq
        r=''.join(reversed(s))
        pos_A = r.find('A')
        pos_T = r.find('T')
        pos_G = r.find('G')
        pos_C = r.find('C')
        if(min(pos_A,pos_T,pos_G,pos_C)>-1 and min(pos_A,pos_T,pos_G,pos_C) < pos_gap):
            pos_gap = len(s) - min(pos_A,pos_T,pos_G,pos_C)

#cut off terminal gaps
    alignment = alignment[:,:pos_gap]
    global max_length
    seq_length = len(alignment[0].seq)
    if max_length<seq_length: max_length = seq_length
    return alignment

#alignment = trim_end_gaps(alignment)


        


# In[5]:


#remove >=50% gap seqs

def rm_gappy(alignment):
    gap_counter = 0.0
    gappy_seqs = []
    ali_len = 0
    f = open('ungapped_alignment','w')
    f.close()
#find sequences with >=50% gaps
    for i in range(len(alignment)):
        for letter in alignment[i].seq:
            if(letter not in['A','T','G','C']):
                gap_counter+=1
        if gap_counter/len(alignment[i].seq)>=0.5:
            gappy_seqs.append(i)
        else:
            SeqIO.write(alignment[i],open('ungapped_alignment','a'),file_format)    
        gap_counter = 0
    alignment = AlignIO.read(open('ungapped_alignment'),file_format)
    return alignment
#alignment = rm_gappy(alignment)





# In[6]:


#reduce to <500 representative seqs
def shrink_alignment(alignment,target_num):
    rep_seqs = 0
    global similarity_threshold
    finished = False
    summary = SummaryInfo(alignment)
    consensus = summary.dumb_consensus(0.5, "-")
    help_shrink(alignment,rep_seqs,similarity_threshold,finished,consensus,target_num)


def help_shrink(alignment, rep_seqs, similarity_threshold, finished, consensus,target_num):
    f = open('representative_seqs','w')
    f.close()
    while not finished:
        for record in alignment:
            hits = 0
            for i in range(len(record.seq)):
                if record.seq[i] == consensus[i]:
                    hits+=1
            if hits/len(record.seq) >= float(similarity_threshold) and rep_seqs<target_num:
                rep_seqs+=1
                SeqIO.write(record,open('representative_seqs','a'),file_format)
        if rep_seqs ==target_num or record.seq == alignment[len(alignment)-1].seq:
            finished=True
            result = alignment
        try:
            alignment = AlignIO.read(open('representative_seqs'),file_format)
            result = alignment
        except:
            similarity_threshold = input('similarity_threshold too high. please type a new one between 0 and 1 that is less than '+str(similarity_threshold)+'\n')
            help_shrink(alignment,0,similarity_threshold,False,consensus,target_num)
    return result
#alignment = shrink_alignment(alignment)


# In[7]:


#make tree
def fasttree(i:int):
    path = file_data[i][0]+'-'+file_data[i][1]+'-'+file_data[i][2]
    cmd = " fasttree -nt "+path+" >" + path+".tree"
    print(subprocess.call(cmd,shell=True,executable='/bin/bash'))

#fasttree(0)


# In[8]:


def fubar(i:int):
    path = file_data[i][0]+'-'+file_data[i][1]+'-'+file_data[i][2]
    cmd = "echo `(echo '1'; echo '4'; echo"+" '"+path+"'; "+"echo "+"'"+path+".tree') | hyphy `"
    print(subprocess.call(cmd,shell=True,executable='/bin/bash'))
#fubar(0)


# In[9]:


def busted(fasta_files,num_seqs):
    f = open('busted_seqs','w')
    f.close()
    global max_length
    for i in fasta_files:
        f1 = open(i,'r')
        f2 = open('busted_seqs','a')
        alignment = AlignIO.read(f1,file_format)
        alignment = trim_end_gaps(alignment)
        alignment = rm_gappy(alignment)
        alignment = shrink_alignment(alignment,num_seqs)
        alignment = AlignIO.read(open('representative_seqs','r'),file_format)

        seq_length = len(alignment[0].seq)
        if max_length>seq_length:
            for j in range(len(alignment)):
                seq = alignment[j].seq
                for i in range(seq_length,max_length):
                    seq = seq+'-'
                    alignment[j].seq = seq
        #print(max_length,seq_length,len(alignment[0].seq))
        AlignIO.write(alignment,f2,file_format)
        f1.close()
        f2.close()
    
    cmd = " fasttree -nt "+path+"/busted_seqs >" + path+"/busted_seqs.tree"
    print(subprocess.call(cmd,shell=True,executable='/bin/bash'))
    
    tree = Tree('busted_seqs.tree',format=1)
    
    found_subtypes = []
    human_clades = []
    clade_ancestors = []
    animal_clade=[]
    #make list of every node of same subtype and find common ancestor of that list. label with {test}. add to list of common ancestors. root on common ancestor of common ancestors.

    #root tree by setting animal lineage as outgroup
    for i in file_data:
        if i[1]!='human':
            animal_lineage=i[2].split('.')[0]

    #populate list of subtypes
    for node in tree.traverse('levelorder'):
        name = node.name.split('.')[0]
        if not name in found_subtypes:
            found_subtypes.append(name)
        if name.casefold()==animal_lineage.casefold():
            animal_clade.append(node)

    #make outgroup = animal clade ancestor
    ancestor = tree.get_common_ancestor(animal_clade)
    tree.set_outgroup(ancestor)
    ancestor.name+='{test}'
    print(tree)

    #list all nodes per human clade
    for i in range(len(found_subtypes)):
        if not found_subtypes[i].casefold()==animal_lineage.casefold():
            clade=[]
            for node in tree.traverse('preorder'):
                name = node.name.split('.')[0]
                if found_subtypes[i].casefold() == name.casefold():
                    clade.append(node)
            human_clades.append(clade)
    #get common ancestor of each clade and label {test}
    for clade in human_clades:
        if len(clade)>1:
            ancestor = tree.get_common_ancestor(clade)
        elif len(clade)==1:
            ancestor = clade[0]
        clade_ancestors.append(ancestor)
        ancestor.name = ancestor.name+'{test}'

    tree.write(outfile='busted_labeled.tree',format=1)
    tree.write(format=1)
    
    path1 = file_data[0][0]+'/busted_seqs'
    path2 = file_data[0][0]+'/busted_labeled.tree'
    cmd = "hyphy busted --alignment busted_seqs --tree busted_labeled.tree --branches test"
    print(subprocess.call(cmd,shell=True,executable='/bin/bash'))
    


# In[10]:


def translate (codon):
    if codon == 'TGG': return 'W'
    elif codon == "TTT" or codon == "TTC": return "F"
    elif codon in ['TTA','TTG','CTT','CTC','CTA','CTG']: return "L"
    elif codon in ['ATT','ATC','ATA']: return "I"
    elif codon == "ATG": return 'M' #START
    elif codon in ['GTT','GTC','GTA','GTG']: return "V"
    elif codon in ['TCT','TCC','TCA','TCG']: return "S"
    elif codon in ['CCT','CCC','CCA','CCG']: return "P"
    elif codon in ['ACT','ACC','ACA','ACG']: return "T"
    elif codon in ['GCT','GCC','GCA','GCG']: return "A"
    elif codon in ['TAT','TAC']: return "Y"
    elif codon in ['TAA','TAG','TGA']: return "-" #STOP
    elif codon in ['CAT','CAC']: return "H"
    elif codon in ['CAA','CAG']: return "Q"
    elif codon in ['AAT','AAC']: return "N"
    elif codon in ['AAA','AAG']: return "K"
    elif codon in ['GAT','GAC']: return "D"
    elif codon in ['GAA','GAG']: return "E"
    elif codon in ['TGT','TGT']: return "C"
    elif codon in ['CGT','CGC','CGA','CGG','AGA','AGG']: return "R"
    elif codon in ['GGT','GGC','GGA','GGG']: return "G"
    else: return '-1'


# In[11]:


#parse json and make array of sites under negative selection
def parse_fubar_json(i:int,consensus):
    path = fasta_files[i]
    fubar_conserved = []
    #summary = SummaryInfo(alignment)
    #consensus = summary.dumb_consensus(0.5, "-")
    with open(fasta_files[i]+'.FUBAR.json','r') as f:
        data = json.load(f)
        for j in data['data partitions']['0']['coverage'][0]:
        #print(i)
        #print(data['MLE']['content']['0'][j]) #site number is in data[data partitions][0][coverage]
            dN = data['MLE']['content']['0'][j][1]
            dS = data['MLE']['content']['0'][j][0]
        #print(dNdS)
            if dN-dS < -0.5:  #header beta-alpha is dN-dS - index 3
            #print(dN-dS,i)
                nt_pos = j*3
                amino = translate(consensus[nt_pos:nt_pos+3])
                fubar_conserved.append([j,amino,dN-dS])
    #print(len(fubar_conserved),fubar_conserved)
    return fubar_conserved

#parse_fubar_json(0)


# In[12]:


file_data = []
max_length=0
for i in range(len(fasta_files)):
    s = [fasta_files[i].split('-')]
    file_data+=s
human_fubars = []
animal_fubars=[]
for i in range(0,len(file_data)):
    gene = file_data[i][0]
    host = file_data[i][1]
    subtype = file_data[i][2].split('.')[0]
    alignment1 = read_alignment(i)
    summary = SummaryInfo(alignment1)
    consensus = summary.dumb_consensus(0.5, "-")    
    alignment1 = trim_end_gaps(alignment1)
    alignment1 = rm_gappy(alignment1)
    alignment1 = shrink_alignment(alignment1,1000)
    fasttree(i)
    fubar(i)
    #time.sleep(20)
    result = [gene[gene.rfind('/')+1:len(gene)],host,subtype,parse_fubar_json(i,consensus)]
    if file_data[i][1]=='human':
        human_fubars+=(result)
    else:
        animal_fubars+=(result)
    


# In[13]:


#find fubar results of same gene diff host where encoded AA is different
diff_conserved = []
animal = animal_fubars[1]
gene = human_fubars[0]
aType = animal_fubars[2]
hType = human_fubars[2]

f = open(aType+'_vs_'+hType+'_differentially_conserved_sites.txt','w')
f.close

if len(human_fubars[3])<len(animal_fubars[3]):
    short=human_fubars[3]
    long = animal_fubars[3]
    s = 'human'
    l = animal
else:
    short = animal_fubars[3]
    long = human_fubars[3]
    s = animal
    l = 'human'

df_s_amino=[]
df_l_amino=[]
df_sites=[]
for i in range(len(short)):
    site = short[i][0]
    amino = short[i][1]
    if amino =='-1': continue
    for j in range(len(long)):
        if long[j][0]==site and long[j][1]!=amino:
            if long[j][1]!='-1':
                result = s+': '+amino+'; '+l+': '+long[j][1]+'; site = '+str(site)+'\n'
                df_s_amino.append(amino)
                df_l_amino.append(long[j][1])
                df_sites.append(site)
                if not '-1' in result:
                    with(open(aType+'_vs_'+hType+'_differentially_conserved_sites.txt','a')) as f:
                        f.write(result)
                    print(result)
data = {s: df_s_amino,
        l: df_l_amino,
        'codon site': df_sites}
df = pd.DataFrame(data)
print(df)


# In[14]:


busted(fasta_files,150) #the int is how many seqs of each lineage you want to analyze with busted


# In[15]:


threshold = 0
diversifying_sites=[]
with open('busted_seqs.BUSTED.json','r') as f:
        data = json.load(f)
        constrained_stat = data['Site Log Likelihood']['constrained'][0]
        for i in range(len(constrained_stat)):
            if float(constrained_stat[i]) >= threshold:
                #print(constrained_stat[i])
                diversifying_sites.append(i)
        df['host adaptive']=['True' if x in diversifying_sites else 'False' for x in df['codon site']]   
        df['constrained statistic'] = [constrained_stat[x] for x in df['codon site']]
with open(gene+'_IHAM_output.csv','w') as f:
    df.to_csv(f)
            


# In[16]:


cmd = 'rm *.json *.tree *.cache *.log ungapped_alignment representative_seqs'
print(subprocess.call(cmd,shell=True,executable='/bin/bash'))
