'''
This script calculates the average number of substititutions at each site in a protein from a BEAST trees file. 
This is done for the entire tree as well as the trunk of the tree.This analysis was done for M1 and NP for 
human and swine influenza. It assumes a 10% burn-in.

Functions
-----------
*``findlengthandterminalnode`` : creates sbatch file
*``gettaxanamedict`` : creates dictionary for replacement taxa names
*``formattree`` : formats tree to newick format so that Bio.Phylo works 
*``readtree`` : uses Bio.Phylo to read newick tree
*``tracetrunk`` : uses Bio.Phylo to trace trunk of last time-stamped taxa to root of tree
*``tracetrunkfromlist`` : performs tracetrunk on a list
*``getmutationsites`` : finds average number of mutations that occur along the tree or trunk
*``writeoutfile`` : writes results 

Input files
-------------
*``prot_aligned.trees`` : BEAST trees

Output files
-------------
*``treeavemutationpersite.csv`` : average mutations at each site for entire tree
*``trunkavemutationpersite.csv`` : average mutations at each site for trunk

'''


import os
import Bio.Phylo
from StringIO import StringIO
import sys
import re

def findlengthandterminalnode(tree):
    '''
    This function returns the number of sequences used to build a beast file *tree* and the length of the protein
    '''
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
    '''
    The beast file *treefile* has the taxa names replaced by a number in the tree - this number taxa name format isn't recognized by biophylo as the parent.
    Thus what this function does is enable one to link the number with a name that is recognized by biophylo. This function creates a dictionary 
    that has the number name as the key and the replacement name as the value. The replacement name is hostnumberinf (host is human and swine).
    This function uses the number of taxa, *terminalnode*, to make replacement names. 
    This function also returns the replacement name of the terminal (last date-stamped) sequence. '''

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
                    count+=1
                    nodedictionary[entry[0]] = '%s%sinf' % (addword,count)                           
    f.close()
    print termnode
    #for nodes in sorted(nodedictionary):
        #print nodes, nodedictionary[nodes]
         
    return nodedictionary, termnode


def formattree(tree,nodedict):
    '''
    This function simplifies and reformats the BEAST trees *tree* so that Bio.Phylo can read them properly. 
    This editing is done assuming a 10% burn-in.

    It removes some information incuding the summary mutation information at the top of the 
    tree and the sequence state information embedded in the tree.

    It does some formatting so that biophlyo can read the tree properly.  
    It takes a treefile that is in the format:
        :[&history_all={31,0.1,K,R}]110.2
    where the [&history...] contains the mutational path and 110.2 is the branch length and swaps the mutation
    history and branchlength so it looks like this:
        :110.2[&history_all={31,0.1,K,R}]
    It replaces the number that represents the taxa with a name from the dictionary *nodedict*, created in gettaxanamedict

    This function returns the edited trees in the BEAST file as a list. 

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

    return last_list

def readtree(tree):
    '''
    This function reads a newick tree *tree* from a string using Bio.Phylo
    '''
    #print 'reading tree'
    readtree = Bio.Phylo.read(StringIO(tree), 'newick')
    return readtree
  
def tracetrunk(tree, terminal):
    '''
    This function reads a newick tree *tree*, and given the name of a terminal node *terminal* returns the path to the root.
    '''
    
    tracepath = Bio.Phylo.BaseTree.TreeMixin.get_path(tree, target = terminal)
    return tracepath

def tracetrunkfromlist(treelist, terminal):
    '''
    This function traces each newick tree in a list of trees *treelist* from a terminal node *terminal* to the root. 
    The result is returned as a list. 
    '''
    tracelist = []
    print "tracing trunk "
    for tree in treelist:
        readt = readtree(tree)
        tracet = tracetrunk(readt, terminal)
        treetrace = [clade.comment for clade in tracet if clade.comment != None]
        treetrace = ', '.join(treetrace)
        tracelist.append(treetrace)
        ###tracelist.append(str(tracet))
    return tracelist
 
def getmutationsites(paths,outfile):
    '''
    This function takes a list of newick trees in *paths* and finds the average number of substitutions that occur at each site.
    This result is returned as a dictionary with aa position as key and average number of substitutions as the value
    '''
    aa_dict = {}
    print "finding aa mutations"
    count = 0
    f=open(outfile,'w')
    for path in paths:
        #print path
        count+=1
        findsubs = re.findall("(\{([\d]*),([\d]*.[\d]*),(\w),(\w)})", path)
        #print findsubs
        for subinfo in findsubs:

            #print subinfo
            #print subinfo[1]
            #print subinfo[2]
            #print subinfo[3]
            #print subinfo[4]
            f.write('%s%s%s:%s,' % (subinfo[3],subinfo[1],subinfo[4],subinfo[2]))
        f.write('\n')

    f.close()
        #print re.findall("(\{([\d]*),([\d]*.[\d]*),(\w),(\w)})", path)
        #for matchobj,aminoacidpos in re.findall("(\{([\d]*),([\d]*.[\d]*),)", path):
            #for matchobj,aminoacidpos in re.findall("(\{([\d]*),([\d]*),(\w),(\w)})", path):
            #print matchobj
           # print aminoacidpos
           # if aminoacidpos not in aa_dict:
             #   aa_dict[aminoacidpos] = 1
            #elif aminoacidpos in aa_dict:
               # aa_dict[aminoacidpos] += 1
    #for aminoacid in aa_dict:
       # aa_dict[aminoacid] = float(aa_dict[aminoacid])/count
       # print aminoacid, aa_dict[aminoacid]
    #return aa_dict
def consensustrunksub(subfile,consensusoutfile):
    #print 'reading from subsummaryfile'
    mastersublist = []
    from collections import defaultdict
    sitedict = defaultdict(list)
    #consensussublist = {}
    f= open(subfile,'r')
    fx = open(consensusoutfile,'w')
    count = 0
    with f as allsubfile:
        for line in allsubfile:
            count+=1
            entry = line.strip().split()
            #print entry
            mastersublist.append(entry)
    f.close()

    cutoff = float(count)*.9
    print "cutoff: %s" % cutoff
    for tree in mastersublist:
        print "sub in tree: %s" % tree
        #checkforduplicationslist = []
        for substitutions in tree:
            print "substitution level %s" % substitutions
            checkforduplicationslist = []
            #print substitution
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
                checkforduplicationslist.append(indivsub[0])
                #print indivsub[0], indivsub[1]
                sitedict[indivsub[0]].append(indivsub[1])
            #print 'length duplist %s' % len(checkforduplicationslist)

            from collections import Counter         
            print [k for k,v in Counter(checkforduplicationslist).items() if v>1]

    for subidentity in sitedict:
        if len(sitedict[subidentity]) > cutoff:
            #print subidentity, sitedict[subidentity]
            blengthtot = 0
            subnumber = 0
            for blength in sitedict[subidentity]:
                blengthtot+=float(blength)
                subnumber+=1
            aveblength = float(blengthtot)/subnumber

            fx.write('%s:%s\n' % (subidentity, aveblength))
        else:
            continue
            #print "not consensus %s %s %s" % (subidentity, sitedict[subidentity], len(sitedict[subidentity]))
            #print subinfo
           # print subinfo[1]
            #print subinfo[1][1]
            #print subinfo[2]
            #sitedict[]
    fx.close()
    


def writeoutfile(aadict, outfile, lengthprot):
    '''
    This function writes the results of the average number of substitutions per site to a csv file *outfile*. 
    *aadict* is a dictionary with amino acid as the key and average number of substitutions as the key. The oufile has 
    a titleline: site, avemutation. Each subsequent line has an amino acid number and the corresponding number of mutations. 
    *lengthprot* is used because if there were no mutations, that site is not in the dictionary, so the value is recorded as
    0 in the outfile. 
    '''

    f = open(outfile, 'w')
    f.write('site,avemutation\n')
    for aa in range(lengthprot):
        aa += 1
        if str(aa) not in aadict:
            f.write('%s,0\n' %(aa))
        else:
            f.write('%s,%s\n' %(aa, aadict[str(aa)]))
    f.close()

def main():
    home = os.path.expanduser("~")

    treefiles = (
                '%s/human/NP/prot_aligned.trees' % os.getcwd(),
                '%s/human/M1/prot_aligned.trees' % os.getcwd(),  
                '%s/swine/NP/prot_aligned.trees' % os.getcwd(),
                '%s/swine/M1/prot_aligned.trees' % os.getcwd(),        
                )

    trunkoutfiles = (
                '%s/human/NP/trunksubidentityandbranchlength.csv'% os.getcwd(),
                '%s/human/M1/trunksubidentityandbranchlength.csv'% os.getcwd(),
                '%s/swine/NP/trunksubidentityandbranchlength.csv'% os.getcwd(),
                '%s/swine/M1/trunksubidentityandbranchlength.csv'% os.getcwd(),
                )
    trunkoutfiles2 = (
                '%s/human/NP/trunksubidentityandlocation.csv'% os.getcwd(),
                '%s/human/M1/trunksubidentityandlocation.csv'% os.getcwd(),
                '%s/swine/NP/trunksubidentityandlocation.csv'% os.getcwd(),
                '%s/swine/M1/trunksubidentityandlocation.csv'% os.getcwd(),
                )

    align_datainput = list(zip(treefiles,trunkoutfiles,trunkoutfiles2)) 
    for entry in align_datainput:
        print entry[0]
        findnodeandlength = True
        if findnodeandlength:
            terminalnode, sequencelength = findlengthandterminalnode(entry[0])
        getnodedict = True
        if getnodedict:
            nodedict,lastnode = gettaxanamedict(entry[0],terminalnode)
        
        treeformat = True
        if treeformat:
            swaptree = formattree(entry[0], nodedict)     

        treemutations = False
        if treemutations:
            treemutations = getmutationsites(swaptree)

        treemutationspersitefile = False
        if treemutationspersitefile:
            writeresults = writeoutfile(treemutations,entry[1], sequencelength)

        tracetrunks = True
        if tracetrunks:
            trunk = tracetrunkfromlist(swaptree, lastnode)

        trunkmutations = True
        if trunkmutations:
            trunkmutations = getmutationsites(trunk,entry[1])

        trunkmutationspersitefile = False
        if trunkmutationspersitefile:
            writeresults = writeoutfile(trunkmutations,entry[1], sequencelength)

        findconsensussub = False
        if findconsensussub:
            consensussub= consensustrunksub(entry[1],entry[2])

        



        


main()