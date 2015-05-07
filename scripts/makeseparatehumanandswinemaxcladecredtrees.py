'''
This script takes a tree that contains both human and swine influenza lineages and creates separate 
human and swine trees. This analysis is done for the maximum clade credibility trees for influenza NP, M1, and PA.

Functions:
------------

Input files:
-------------
*``maxcladecredibility.tre`` : maximum clade credibility tree that contains human and swine sequences

Output files:
--------------
*``humanmaxcladecredibility.tre`` : maximum clade credibility tree with only human sequences
*``swinemaxcladecredibility.tre`` : maximum clade credibility tree with only swine sequences

'''
#steps:
#read in appropriate tree
#get human list and swine list - from file read-in
#reformat tree
#prune trees

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

    f = open(treefile, 'r')
    with f as treef:
        for line in treef:
            entry = line.strip().split()
            if len(entry) == 2:
                if (entry[0]) in nodelist: 
                    entry[1] = entry[1].replace("'", "")
                    entry[1] = entry[1].replace(",", "")
                    if ('swine' or 'Swine') in entry[1]:
                        swine +=1
                        nodedictionary[entry[0]] = 'swine%sinf' % swine
                    else:
                        human +=1
                        nodedictionary[entry[0]] = 'human%sinf' % human              
    f.close()

    humanterminalnode = ''
    swineterminalnode = ''
    j = nodelistlength
    for i in range(nodelistlength):
        if humanterminalnode =='':
            if 'swine' not in nodedictionary[str(j-i)]:
                humanterminalnode += nodedictionary[str(j-i)]
        if swineterminalnode =='':
            if 'swine' in nodedictionary[str(j-i)]:
                swineterminalnode += nodedictionary[str(j-i)]
        if swineterminalnode != '' and humanterminalnode != '':
            break
       
    print human, swine, humanterminalnode, swineterminalnode
         
    return nodedictionary,humanterminalnode,swineterminalnode

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

def main():
    treefiles = (
        '%s/human/M1/prot_aligned.trees' % os.getcwd(),  
        '%s/human/NP/prot_aligned.trees' % os.getcwd(),  
        '%s/human/PA/prot_aligned.trees' % os.getcwd(),
        )
    maxcladecredibilitytreenumber = (21310000,7610000,8060000)

    humantreeoutfiles = (
        '%s/human/M1/prot_aligned.trees' % os.getcwd(),  
        '%s/human/NP/prot_aligned.trees' % os.getcwd(),  
        '%s/human/PA/prot_aligned.trees' % os.getcwd(),
        )
    
    align_datainput = list(zip(treefiles,humantreeoutfiles,swinetreeoutfiles,humantrunkoutfiles,swinetrunkoutfiles)) #terminalnodenumberlabel,terminalnodereplacementlabel, outfiles, len_protein))



main()