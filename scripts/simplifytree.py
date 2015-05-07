'''
1. The BEAST maxcladecredibility tree does not work as is for datamonkeyan online server for dnds analysis. Simplify the BEAST maxcladecredibility tree using ``simplifytree.py``. 

Input files 
------------
*``maxcladecredibility.tre`` : maximum clade credibility tree
*``cds_aligned.fasta`` : aligned cDNA sequences used for making tree file

Output files
--------------
*``testsimplenexus.tre`` : output nexus tree for datamonkey
'''

import os
import Bio.Phylo
from StringIO import StringIO
import sys
import re
import listFASTA
from Bio import SeqIO
from collections import OrderedDict

def findlengthandterminalnode(tree):
    terminalnode = ''
    lengthseq = 0
    f = open(tree, 'r')
    with f as treef:
        for line in treef:  
            if terminalnode == '':
                matchnode = re.findall('Dimensions ntax=(\d+)',line)
                if matchnode:
                    terminalnode += matchnode[0]
                    print "node %s" % terminalnode
            if lengthseq == 0:
                findlength = re.findall('\[&states="([A-Z]+)"\]', line)
                if findlength:
                    lengthseq += len(str(findlength[0]))
                    print "lengthseq %s" % lengthseq
            else:
                f.close()
                break
    return terminalnode, lengthseq



def gettaxanamedict(treefile,terminalnode):
    '''the beast file has the taxa names replaced by a number - this number taxa name format isn't recognized by biophylo as the parent.
    Thus what this function does is enable one to replace the number with the actual taxa name for all of the beast trees. This function retrieves 
    the number and the corresponding taxa and puts it in a dictionary with the number as a key and the taxa name as the value.'''

    nodedictionary = {}
    nodelist = [str(i+1) for i in range(int(terminalnode))]
    nodelistlength = len(nodelist)
    taxaabreviationdict = OrderedDict()
    human = 0
    swine = 0
    count = 0
    if 'swine' in treefile:
        addword ='swine'
        termnode = 'swine%sinf' % terminalnode
    else:
        addword ='human'
        termnode = 'human%sinf' % terminalnode


    f = open(treefile, 'r')
    with f as treef:
        for line in treef:
            entry = line.strip().split()
            if len(entry) == 2:
                if (entry[0]) in nodelist: 
                    entry[1] = entry[1].replace("'", "")
                    entry[1] = entry[1].replace(",", "")
                    #if ('swine' or 'Swine' or 'SW' or 'sw') in entry[1]:
                    count+=1
                    #swine +=1
                    nodedictionary[entry[0]] = '%s%sinf' % (addword,count)
                    taxaabreviationdict[entry[1]] = '%s%sinf' % (addword,count)
                   # else:
                       # human +=1
                       # nodedictionary[entry[0]] = 'human%sinf' % human              
    f.close()

   # if 'swine' in treefile:
      #  termnode = 'swine%sinf' % terminalnode
   # else:
       # termnode = 'human%sinf' % terminalnode

    print termnode
    for nodes in sorted(nodedictionary):
        print nodes, nodedictionary[nodes]
         
    return nodedictionary, termnode, taxaabreviationdict


def swapbranchlengthandhistory(tree,nodedict):
    '''
    this function takes a treefile that is in the format:

    :[&history_all={31,0.1,K,R}]110.2

    where the [&history...] contains the mutational path and 110.2 is the branch length

    the output of this function swaps the mutation history and branchlength so it looks like this:

    :110.2[&history_all={31,0.1,K,R}]

    this format is then used for reading the tree using biopython

    '''
    print "reading file and formatting"
    list_trees = []
    list_remove = []
    new_list = []
    list_a = []
    if not os.path.isfile(tree):
        raise IOError('Failed to find infile of %s' % tree)
    f = open(tree, 'r')
    with f as treef:
        for line in treef:
            match = re.findall("tree TREE1", line)
            for state in match:
               # if (int(state) > 3000000) :
                no_top = re.sub("tree TREE1 = \[&R\] ","", line)
                
                list_a.append(no_top)
    f.close()
    
    for tree in list_a:
        no_states = re.sub('(\[&states.set={"[A-Z,""]*"}[\%_\w\d,"&-=\.{}]*\])', "", tree)
        list_trees.append(no_states)
   # print list_trees
  
    for tree in list_trees:
        for matchobj in re.finditer("(\[&history_all=\{[{}\w,-\.]*\}\])(\d+\.?\d*)", tree):
            tree = tree[:matchobj.start(0)] + tree[matchobj.end(1):matchobj.end(0)] + tree[matchobj.start(0):matchobj.end(1)] + tree[matchobj.end(0):]
        new_list.append(tree)
  
    last_list = []
    i = 0
    for tree in new_list:
        i +=1
        if i%100 == 0:
            print i
        matchobj = re.findall("((,|\()(\d+):)", tree)
        for match in matchobj:
            #print nodedict
          #  print match[2]
           # print nodedict[1]
           # print nodedict[str(match[2])]
            replacement = "%s:" % (nodedict[match[2]])
            tree = re.sub('(?<=[,\(])%s:'%match[2],replacement, tree)
        last_list.append(tree)
    print last_list
    for item in last_list:
        item = item.replace('\n','')

    print item

    return item #last_list

def readtree(tree):
    '''
    This function reads a newick tree from a string 
    '''
    #print 'reading tree'
    readtree = Bio.Phylo.read(StringIO(tree), 'newick')
    print readtree
    return readtree

def testtree(treefilename):


    f = open(treefilename, 'r')
    with f as treef:
        for line in treef:
            entry = line.strip().split()
            tree = entry
            print entry
    readtree = Bio.Phylo.read(StringIO(tree), 'newick')
    print readtree

def sequenceread(outfile,seqfile,abbreviationdict,terminalnode,tree):
    sequences = listFASTA.listFASTA(seqfile)
    f = open(outfile,'w')
    f.write('#NEXUS\n')
    f.write('BEGIN TAXA;\n')
    f.write('   DIMENSIONS NTAX = %s;\n'% len(sequences))
    f.write('   TAXLABELS\n')
    print terminalnode
    for entry in abbreviationdict:
        if abbreviationdict[entry] ==terminalnode:
            f.write("       '%s' ;\n" %abbreviationdict[entry])
        else: # abbreviationdict[entry] != terminalnode:
            f.write("       '%s'\n" %abbreviationdict[entry])
        #if abbreviationdict[entry] =='human207inf'
           # f.write("       '%s' ;\n" %abbreviationdict[entry])
    f.write('END;\n')
    f.write('BEGIN CHARACTERS;\n')
    f.write('   DIMENSIONS NCHAR = %s;\n' % len(sequences[0].seq))
    f.write('   FORMAT\n')
    f.write('       DATATYPE = DNA\n')
    f.write('       GAP=-\n')
    f.write('       MISSING=?\n')
    f.write('   ;\n')
    f.write('MATRIX\n')
    orderedseq = []
    numberseq = len(sequences)
    count =0
    for i,entry in enumerate(abbreviationdict): #seqs in sequences:
        for seqs in sequences:
            if entry == seqs.description:
                count +=1
                if i != numberseq -1:
                    #count +=1
               # f
                    f.write("   '%s'     %s\n" %(abbreviationdict[entry], seqs.seq))
                else:
                    f.write("   '%s'     %s;\n" %(abbreviationdict[entry], seqs.seq))
    f.write('END;\n')
    f.write('BEGIN TREES;\n')
    f.write('   TREE tree = %s\n' %(tree))
    f.write('END;\n')
    print count
    f.close()

def main():
    home = os.path.expanduser("~")

    testread = testtree('%s/testtree.tre'% os.getcwd())
    
    treefiles = (
               # '%s/human/NP/maxcladecredibility.tre' % os.getcwd(), 
               # '%s/human/M1/maxcladecredibility.tre' % os.getcwd(),  
               # '%s/swine/NP/maxcladecredibility.tre' % os.getcwd(),
               # '%s/swine/M1/maxcladecredibility.tre' % os.getcwd(),  
                '%s/human/HA_H3/maxcladecredibility.tre' % os.getcwd(), 
                '%s/swine/HA_H3/maxcladecredibility.tre' % os.getcwd(),         
                )

    treeoutfiles = (
               # '%s/human/NP/testsimplenexus.tre'% os.getcwd(),
               # '%s/human/M1/testsimplenexus.tre'% os.getcwd(),
               # '%s/swine/NP/testsimplenexus.tre'% os.getcwd(),
               # '%s/swine/M1/testsimplenexus.tre'% os.getcwd(), 
                '%s/human/HA_H3/testsimplenexus.tre'% os.getcwd(),
                '%s/swine/HA_H3/testsimplenexus.tre'% os.getcwd(),          
                )

    alignedseqs = (
               # '%s/human/NP/cds_aligned.fasta'% os.getcwd(),
               # '%s/human/M1/cds_aligned.fasta'% os.getcwd(),
               # '%s/swine/NP/cds_aligned.fasta'% os.getcwd(),
               # '%s/swine/M1/cds_aligned.fasta'% os.getcwd(),
                '%s/human/HA_H3/cds_aligned.fasta'% os.getcwd(),
                '%s/swine/HA_H3/cds_aligned.fasta'% os.getcwd(),
                )


    align_datainput = list(zip(treefiles,treeoutfiles,alignedseqs)) #terminalnodenumberlabel,terminalnodereplacementlabel, outfiles, len_protein))
    for entry in align_datainput:
        print entry[0]
        findnodeandlength = True
        if findnodeandlength:
            terminalnode, sequencelength = findlengthandterminalnode(entry[0])
        getnodedict = True
        if getnodedict:
            nodedict,lastnode,abbreviationdict = gettaxanamedict(entry[0],terminalnode)
       # readseq = True
       # if readseq:
          #  seqdict = sequenceread()
        
        swapbranchandhistory = True
        if swapbranchandhistory:
            swaptree = swapbranchlengthandhistory(entry[0], nodedict) 

        output = True
        
        if output:

            out = sequenceread(entry[1],entry[2],abbreviationdict,lastnode,swaptree)

        readsimpletree = False
        if readsimpletree:
            importedtree = readtree(swaptree)
    
main()







