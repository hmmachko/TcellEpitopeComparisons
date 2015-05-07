'''
This script calculates the 
'''
#this script tries to read along the trunk of a nexus tree from the most recent taxa to the root
#and extract the mutations that occur along the trunk. 

import os
import Bio.Phylo
from StringIO import StringIO
import sys
import re

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
         
    return nodedictionary, termnode


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
            match = re.findall("tree STATE_(\d+)", line)
            for state in match:
                if (int(state) > 3000000) :
                    no_top = re.sub("tree STATE_[^()]*","", line)
                    list_a.append(no_top)
    f.close()
    
    for tree in list_a:
        no_states = re.sub('(\[&states="[A-Z]+"\])', "", tree)
        list_trees.append(no_states)
  
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
            replacement = "%s:" % (nodedict[match[2]])
            tree = re.sub('(?<=[,\(])%s:'%match[2],replacement, tree)
        last_list.append(tree)
    #print last_list

    return last_list

def readtree(tree):
    '''
    This function reads a newick tree from a string 
    '''
    #print 'reading tree'
    readtree = Bio.Phylo.read(StringIO(tree), 'newick')
    return readtree
  
def tracetrunk(tree, terminal):
    
    #This function reads a newick tree and given the name of a terminal node returns the path to the root, does not include root
    
    tracepath = Bio.Phylo.BaseTree.TreeMixin.get_path(tree, target = terminal)
    #print tracepath
    return tracepath

def tracetrunkfromlist(treelist, terminal):

    tracelist = []
    print "tracing trunk "
    for tree in treelist:
        readt = readtree(tree)
        #print readt
        #print terminal
       # print "trace of trunk"
        tracet = tracetrunk(readt, terminal)
       # print tracet
        tracelist.append(str(tracet))

    return tracelist
 
def getmutationsites(paths):
    '''
    This function returns dictionary with aa position as key and number of mutations as value
    '''
    aa_dict = {}
    print "finding aa mutations"
    count = 0
    for path in paths:
        count+=1
        for matchobj,aminoacidpos in re.findall("(\{([\d]*),)", path):
            if aminoacidpos not in aa_dict:
                aa_dict[aminoacidpos] = 1
            elif aminoacidpos in aa_dict:
                aa_dict[aminoacidpos] += 1
    for aminoacid in aa_dict:
        aa_dict[aminoacid] = float(aa_dict[aminoacid])/count
        print aminoacid, aa_dict[aminoacid]
    return aa_dict

def writeoutfile(aadict, outfile, lengthprot):

    f = open(outfile, 'w')
    f.write('site,avemutation\n')
    for aa in range(lengthprot):
        aa += 1
        if str(aa) not in aadict:
            f.write('%s,0\n' %(aa))
        else:
            f.write('%s,%s\n' %(aa, aadict[str(aa)]))
    f.close()

def MakeSubTree(treelist, prune):
    '''This function takes the tree from maintreefile, and prunes
    any terminal nodes that appear in the list prune. The new subtree
    is then saved as a newick tree file. The function uses methods from
    the Bio.Phylo module.

    treefile -> File name for tree written in Newick format. 
    subtreefile -> File name for tree generated in this function. This tree
       will have the nodes in prune pruned.
    subtreealignmentfile -> File name for alignment of sequences 
       in the subtree.
    seqs -> A dictionary where the key is the sequence name and the value
       is the sequence. This contains the aligned sequences used to make
       the tree from treefile.
    prune -> A list of the names of terminal nodes in the tree 
       that will get pruned.
    '''

    print 'pruning'
    pruned_trees = []
    i =0
    for entry in treelist:
        i +=1
        if i%100 == 0:
            print i
        tree = readtree(entry)
        term = tree.get_terminals()
        #print tree

        for nodename in prune:
            node=[node for node in tree.find_clades(name='.*%s.*' % nodename)]
            #print nodename, 
            #print node,node[0]
            tree.prune(node[0])
        treeformat = Bio.Phylo.BaseTree.Tree.format(tree,'newick')
        pruned_trees.append(treeformat)
    print len(pruned_trees)

    return pruned_trees


def PruneList(nodedict):
    '''This function returns a list of node names to prune from the tree. The
    list of nodes is the complement to the list of nodes whose names match
    user-specified patterns. For example, if the pattern is 'H3N2', the list
    of nodes to prune will be nodes from H1N1 and H2N2 viral subtypes.

    treefile -> File name for tree written in Newick format. 
    pattern -> A list of substrings that are found in the names of the terminal
       nodes in a phylogenetic tree. For example, one element of pattern could
       be H1N1, and this would match terminal nodes that had names like 
       1984.62_36_A/Managua/3759.02/2008_HOST_Human_H1N1. Pattern could be set to
       a set of years (1918-2014), or viral subtypes (H1N1, H2N2, H3N2).
    problem_symbol -> There is a terminal node name, 957.50_132_A/Leningrad/134/47/ts+18/1957_HOST_Human_H2N2,
       which has a + in it. This interferes with looking for node names in the tree, so the
       PruneList function has code to deal with this specific name.
    prune -> A list of node names to prune from the tree. This list is returned by the function.
    '''

    prune_swine = []
    prune_human = []
    
    for key in nodedict:
        print nodedict[key]
        if 'swine' in nodedict[key]:
            prune_swine.append(nodedict[key])
        else:
            prune_human.append(nodedict[key])

    return prune_swine,prune_human



def main():
    home = os.path.expanduser("~")

    treefiles = (
                #'%s/human/NP/prot_aligned.trees' % os.getcwd(),
                #'%s/human/M1/prot_aligned.trees' % os.getcwd(),  
                #'%s/swine/NP/prot_aligned.trees' % os.getcwd(),
                #'%s/swine/M1/prot_aligned.trees' % os.getcwd(), 
                '%s/human/HA_H3/prot_aligned.trees' % os.getcwd(),
                '%s/swine/HA_H3/prot_aligned.trees' % os.getcwd(),           
                )

    treeoutfiles = (
                #'%s/human/NP/treeavemutationpersite.csv'% os.getcwd(),
               # '%s/human/M1/treeavemutationpersite.csv'% os.getcwd(),
               # '%s/swine/NP/treeavemutationpersite.csv'% os.getcwd(),
               # '%s/swine/M1/treeavemutationpersite.csv'% os.getcwd(), 
                '%s/human/HA_H3/treeavemutationpersite.csv'% os.getcwd(),
                '%s/swine/HA_H3/treeavemutationpersite.csv'% os.getcwd(),          
                )

    trunkoutfiles = (
               # '%s/human/NP/trunkavemutationpersite.csv'% os.getcwd(),
               # '%s/human/M1/trunkavemutationpersite.csv'% os.getcwd(),
               # '%s/swine/NP/trunkavemutationpersite.csv'% os.getcwd(),
               # '%s/swine/M1/trunkavemutationpersite.csv'% os.getcwd(),
                '%s/human/HA_H3/trunkavemutationpersite.csv'% os.getcwd(),
                '%s/swine/HA_H3/trunkavemutationpersite.csv'% os.getcwd(),
                )


    align_datainput = list(zip(treefiles,treeoutfiles,trunkoutfiles)) #terminalnodenumberlabel,terminalnodereplacementlabel, outfiles, len_protein))
    for entry in align_datainput:
        print entry[0]
        findnodeandlength = True
        if findnodeandlength:
            terminalnode, sequencelength = findlengthandterminalnode(entry[0])
        getnodedict = True
        if getnodedict:
            nodedict,lastnode = gettaxanamedict(entry[0],terminalnode)
        
        swapbranchandhistory = True
        if swapbranchandhistory:
            swaptree = swapbranchlengthandhistory(entry[0], nodedict)     

        getprunelist = False
        if getprunelist:
            human, swine = PruneList(nodedict)
        getsubtree = False
        if getsubtree:
            subtree = MakeSubTree(swap, human)
           # swinesubtree = MakeSubTree(swap, swine)

        treemutations = True
        if treemutations:
            treemutations = getmutationsites(swaptree)
            #swinetreemutations = getmutationsites(swinesubtree)

        treemutationspersitefile = True
        if treemutationspersitefile:
            writeresults = writeoutfile(treemutations,entry[1], sequencelength)
            #writeresults = writeoutfile(swinetreemutations,entry[2], sequencelength)

        tracetrunks = True
        if tracetrunks:
            trunk = tracetrunkfromlist(swaptree, lastnode)
            #swinetrunk = tracetrunkfromlist(swinesubtree, swinetermnode)

        trunkmutations = True
        if trunkmutations:
            trunkmutations = getmutationsites(trunk)
          #  swinetrunkmutations = getmutationsites(swinetrunk)

        trunkmutationspersitefile = True
        if trunkmutationspersitefile:
            writeresults = writeoutfile(trunkmutations,entry[2], sequencelength)
           # writeresults = writeoutfile(swinetrunkmutations,entry[4], sequencelength)

        testtrace = False
        if testtrace:
            trace = tracetrunkfromlist(swap,humantermnode)


        


main()