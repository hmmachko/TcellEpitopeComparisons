'''
This script summarizes the substitutions that occur along the trunk of at least 90% of the BEAST trees (excluding the 10% burnin), and the year that the substition occured. 
If there is a substitution that happens twice in a few trees, we take the substitution branch length that is closest to the average branch length 
calculated from all of the trees that only had a single occurance of that substitution. There is one substitution that occurs twice:
swine NP S473N (which occur at very different times), and is separated by if it occurs above or below a branch length of 30.   
This also converts branch length to year. The last year is 2013 for all lineages. The branch length is subtracted from the year and the resulting
year is taken (year is rounded down).

Functions
-----------
*``findlengthandterminalnode`` : creates sbatch file
*``gettaxanamedict`` : creates dictionary for replacement taxa names
*``formattree`` : formats tree to newick format so that Bio.Phylo works 
*``tracetrunk`` : uses Bio.Phylo to trace trunk of last time-stamped taxa to root of tree
*``tracetrunkfromlist`` : performs tracetrunk on a list
*``getmutationsites`` : finds average number of mutations that occur along the tree or trunk
*``consensustrunksub`` : records substitutions that occur along the trunk of at least 90% of the trees
*``consensustrunksubswinenp`` : records substitutions that occur along the trunk of at least 90% of the trees for swine NP, which has a substitution that occurs twice.
*``convertbranchlengthtoyear`` : converts the branch length to the year and records substitution and year for supplemental file

Input files
-------------
*``prot_aligned.trees`` : BEAST trees

Output files
-------------
*``treeavemutationpersite.csv`` : average mutations at each site for entire tree
*``trunkavemutationpersite.csv`` : average mutations at each site for trunk
*``trunksubidentityandbranchlength.csv``: all trunk substitutions and corresponding branch length for each tree
*``trunksubstitutionsummary.csv`` : trunk substitutions that occur in at least 90% of trees and branch length
*``trunksubstitutionyearsummary.csv`` : trunk substitions that occur in at least 90% of trees and year

'''


import os
import Bio.Phylo
from StringIO import StringIO
import sys
import re
import numpy as np

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

            f.write('%s%s%s:%s,' % (subinfo[3],subinfo[1],subinfo[4],subinfo[2]))
        f.write('\n')

    f.close()

def consensustrunksub(subfile,consensusoutfile):
    '''
    this function takes a *subfile*, which is a file where each line contains all of the substitutions that occur along the trunk of
    a BEAST tree. The format is like this: A45G:25.4 (where the number after the colon is the branchlength). 
    This function then determines which substitutions occur in at least 90% of trees, and records them in *conesnsusoutfile* in the same format
    as above. 
    ''' 
    #print 'reading from subsummaryfile'
    mastersublist = []
    duplicationcounts = {}
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

        for substitutions in tree:
            #print "substitution level %s" % substitutions
            duplicationsinsingletree = []
            checkforduplicationslist = []
            #print substitution
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
                checkforduplicationslist.append(indivsub[0])
            from collections import Counter 
            dup =  [k for k,v in Counter(checkforduplicationslist).items() if v>1]       
            if len(dup) >0:
                #print dup
                for mut in dup:
                    if mut not in duplicationcounts:
                        duplicationcounts[mut] = 1
                    else:
                        duplicationcounts[mut]+=1
            for indivsub in subinfo:
                if indivsub[0] not in dup:
                    sitedict[indivsub[0]].append(indivsub[1]) ###move down in script
    for duplication in duplicationcounts:
        print duplication, duplicationcounts[duplication]

    initialsiteaverageblength = {}
    #find site average branch length discluding any duplications
    for subidentity in sitedict:
        blengthtot = 0
        subnumber = 0
        for blength in sitedict[subidentity]:
            blengthtot+=float(blength)
            subnumber+=1
        aveblength = float(blengthtot)/subnumber
        initialsiteaverageblength[subidentity] = aveblength

    #now go back and only include the duplication that is closest to the initial average
    print 'initial averages'
    #print initialsiteaverageblength
    for tree in mastersublist:
        #print "sub in tree: %s" % tree
        #checkforduplicationslist = []
        for substitutions in tree:
            #print "substitution level %s" % substitutions
            checkforduplicationslist2 = []
            singletreebranchlengths = []
            #print substitution
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
                checkforduplicationslist2.append(indivsub[0])
                singletreebranchlengths.append(indivsub[1])

            D = defaultdict(list)
            for i,item in enumerate(checkforduplicationslist2):
                D[item].append(i)
            D = {k:v for k,v in D.items() if len(v)>1}  
            #print D 
            #print singletreebranchlengths    
            if len(D) >0:
                #print dup
                for mut in D:
                    distancefromaverage = []
                    #print mut
                    for branchl_index in D[mut]:
                        #print 'entry of dup lengths %s ' % D[mut]
                    
                        #print 'average %s ' % initialsiteaverageblength[mut]
                        #print 'dup length %s ' % singletreebranchlengths[branchl_index]
                        distfromave = float(initialsiteaverageblength[mut]) - float(singletreebranchlengths[branchl_index])
                        distancefromaverage.append(abs(distfromave))
                        #print 'dist from ave %s ' % distfromave
                    min_dist = np.min(distancefromaverage)
                   #print min_dist
                    #print index_for_min_dist
                    #print sitedict[mut]
                    #print D[mut][index_for_min_dist]
                    #print singletreebranchlengths[D[mut][index_for_min_dist]]
                    sitedict[mut].append(min_dist)
                    #sitedict[mut].append(singletreebranchlengths[D[mut][index_for_min_dist]])
                    #print sitedict[mut]
    finaldict = {}                
    for subidentity in sitedict:
        if len(sitedict[subidentity]) > cutoff:
            #print subidentity, sitedict[subidentity]
            blengthtot = 0
            subnumber = 0
            for blength in sitedict[subidentity]:
                blengthtot+=float(blength)
                subnumber+=1
            aveblength = float(blengthtot)/subnumber
            finaldict[subidentity] = aveblength


            #fx.write('%s:%s\n' % (subidentity, aveblength))
        else:
            continue

    for w in sorted(finaldict, key=finaldict.get, reverse=True):
        fx.write('%s:%s\n' % (w, finaldict[w]))
    fx.close()
    


def consensustrunksubswinenp(subfile,consensusoutfile):
    #print 'reading from subsummaryfile'
    mastersublist = []
    duplicationcounts = {}
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

    S473N_high = []
    S473N_low = []

    cutoff = float(count)*.9
    print "cutoff: %s" % cutoff
    for tree in mastersublist:

        for substitutions in tree:
            #print "substitution level %s" % substitutions
            duplicationsinsingletree = []
            checkforduplicationslist = []
            #print substitution
            counthigh=0
            countlow=0
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
   
                if indivsub[0] == 'S473N':
                    if float(indivsub[1]) > 30:
                        counthigh+=1
                        S473N_high.append(float(indivsub[1]))
                    else:
                        countlow+=1

                        S473N_low.append(float(indivsub[1]))
            if countlow >=2:
                print 'more than 2 low'
            if counthigh >=2:
                print 'more than 2 high' 
    #print S473N_high
    print len(S473N_high)
    print len(S473N_low)
    print np.min(S473N_high) 
            
    print np.max(S473N_high)
    print np.min(S473N_low) 
    print np.max(S473N_low)
    print 'means'
    print np.mean(S473N_high)
    print np.mean(S473N_low)

    #fx.write('S473N:%s\n' % (np.mean(S473N_high)))
    #fx.write('S473N:%s\n' % (np.mean(S473N_low)))
                    #print indivsub[1]


    for tree in mastersublist:

        for substitutions in tree:
            #print "substitution level %s" % substitutions
            duplicationsinsingletree = []
            checkforduplicationslist = []
            #print substitution
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
                checkforduplicationslist.append(indivsub[0])
            from collections import Counter 
            dup =  [k for k,v in Counter(checkforduplicationslist).items() if v>1]       
            if len(dup) >0:
                #print dup
                for mut in dup:
                    if mut not in duplicationcounts:
                        duplicationcounts[mut] = 1
                    else:
                        duplicationcounts[mut]+=1
            for indivsub in subinfo:
                if indivsub[0] not in dup:
                    sitedict[indivsub[0]].append(indivsub[1]) ###move down in script
    for duplication in duplicationcounts:
        print duplication, duplicationcounts[duplication]

    initialsiteaverageblength = {}
    #find site average branch length discluding any duplications
    for subidentity in sitedict:
        blengthtot = 0
        subnumber = 0
        for blength in sitedict[subidentity]:
            blengthtot+=float(blength)
            subnumber+=1
        aveblength = float(blengthtot)/subnumber
        initialsiteaverageblength[subidentity] = aveblength

    #now go back and only include the duplication that is closest to the initial average
    print 'initial averages'
   # print initialsiteaverageblength
    for tree in mastersublist:
        #print "sub in tree: %s" % tree
        #checkforduplicationslist = []
        for substitutions in tree:
            #print "substitution level %s" % substitutions
            checkforduplicationslist2 = []
            singletreebranchlengths = []
            #print substitution
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", substitutions)
            for indivsub in subinfo:
                checkforduplicationslist2.append(indivsub[0])
                singletreebranchlengths.append(indivsub[1])

            D = defaultdict(list)
            for i,item in enumerate(checkforduplicationslist2):
                D[item].append(i)
            D = {k:v for k,v in D.items() if len(v)>1}  
            #print D 
            #print singletreebranchlengths    
            if len(D) >0:
                #print dup
                for mut in D:
                    distancefromaverage = []
                    #print mut
                    for branchl_index in D[mut]:
                        #print 'entry of dup lengths %s ' % D[mut]
                    
                        #print 'average %s ' % initialsiteaverageblength[mut]
                        #print 'dup length %s ' % singletreebranchlengths[branchl_index]
                        distfromave = float(initialsiteaverageblength[mut]) - float(singletreebranchlengths[branchl_index])
                        distancefromaverage.append(abs(distfromave))
                        #print 'dist from ave %s ' % distfromave
                    min_dist = np.min(distancefromaverage)
                    #print min_dist
                    #print index_for_min_dist
                    #print sitedict[mut]
                    #print D[mut][index_for_min_dist]
                    #print singletreebranchlengths[D[mut][index_for_min_dist]]
                    sitedict[mut].append(min_dist)
                    #sitedict[mut].append(singletreebranchlengths[D[mut][index_for_min_dist]])
                    #print sitedict[mut]
    finaldict = {}
    for subidentity in sitedict:
        if len(sitedict[subidentity]) > cutoff:
            #print subidentity, sitedict[subidentity]
            blengthtot = 0
            subnumber = 0
            for blength in sitedict[subidentity]:
                blengthtot+=float(blength)
                subnumber+=1
            aveblength = float(blengthtot)/subnumber

            if subidentity != 'S473N':
                finaldict[subidentity] = aveblength


                #fx.write('%s:%s\n' % (subidentity, aveblength))
        else:
            continue
    finaldict['S473N_1'] = np.mean(S473N_high)
    finaldict['S473N_2'] = np.mean(S473N_low)
    #fx.write('S473N:%s\n' % (np.mean(S473N_high)))
    #fx.write('S473N:%s\n' % (np.mean(S473N_low)))
    for w in sorted(finaldict, key=finaldict.get, reverse=True):
        if w == 'S473N_1':
            fx.write('S473N:%s\n' % (finaldict[w]))
        elif w == 'S473N_2':
            fx.write('S473N:%s\n' % (finaldict[w]))
        else:
            fx.write('%s:%s\n' % (w, finaldict[w]))
        #print w, d[w]
    fx.close()

def convertbranchlengthtoyear(infile,outfile):
    
    mostrecent = 2013
    listsubandyear = []

    f= open(infile, 'r')
    fx=open(outfile, 'w')
    with f as subfile:
        for line in subfile:
            entry = line.strip().split()
            entry = str(entry)
            subinfo = re.findall("([\w\d]*):([\d]*.[\d]*)", entry)
            for indivsub in subinfo:
                year = mostrecent - float(indivsub[1])
                fx.write('%s,%d\n' % (indivsub[0], int(year)))
                listsubandyear.append([indivsub[0],year])


    f.close()
    fx.close()
    return listsubandyear

def combinefiles(humanM1, humanNP, swineM1,swineNP, outfile):
    f= open(outfile, 'w')
    f.write('summary of substitutions that occur along the trunk of at least 90% of trees\n\n')
    f.write('human influenza M1 substitution,year,,human influenza NP substitution,year,,swine influenza M1 substitution,year,,swine influenza NP substitution,year\n')
    for i, entry in enumerate(humanNP):
        print i,entry
        if i <=(len(humanM1)-1):
            f.write('%s,%s,,' % (humanM1[i][0],int(humanM1[i][1])))
        else:
            f.write(',,,')
        f.write('%s,%s,,' % (humanNP[i][0],int(humanNP[i][1])))
        if i <=(len(swineM1)-1):
            f.write('%s,%s,,' % (swineM1[i][0],int(swineM1[i][1])))
        else:
            f.write(',,,')
        if i <=(len(swineNP)-1):
            f.write('%s,%s,,\n' % (swineNP[i][0],int(swineNP[i][1])))
        else:
            f.write(',,,\n')




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
                '%s/human/NP/trunksubstitutionsummary.csv'% os.getcwd(),
                '%s/human/M1/trunksubstitutionsummary.csv'% os.getcwd(),
                '%s/swine/NP/trunksubstitutionsummary.csv'% os.getcwd(),
                '%s/swine/M1/trunksubstitutionsummary.csv'% os.getcwd(),
                )
    trunkoutfileswithyear = (
                '%s/human/NP/trunksubstitutionyearsummary.csv'% os.getcwd(),
                '%s/human/M1/trunksubstitutionyearsummary.csv'% os.getcwd(),
                '%s/swine/NP/trunksubstitutionyearsummary.csv'% os.getcwd(),
                '%s/swine/M1/trunksubstitutionyearsummary.csv'% os.getcwd(),
        )

    align_datainput = list(zip(treefiles,trunkoutfiles,trunkoutfiles2,trunkoutfileswithyear)) 
    for entry in align_datainput:
        print entry[0]
        findnodeandlength = False
        if findnodeandlength:
            terminalnode, sequencelength = findlengthandterminalnode(entry[0])
        getnodedict = False
        if getnodedict:
            nodedict,lastnode = gettaxanamedict(entry[0],terminalnode)
        
        treeformat = False
        if treeformat:
            swaptree = formattree(entry[0], nodedict)     

        tracetrunks = False
        if tracetrunks:
            trunk = tracetrunkfromlist(swaptree, lastnode)

        trunkmutations = False
        if trunkmutations:
            trunkmutations = getmutationsites(trunk,entry[1])


        findconsensussub = False
        if findconsensussub:
            if entry[0] == '%s/swine/NP/prot_aligned.trees' % os.getcwd():
                consensusswinenp = consensustrunksubswinenp(entry[1],entry[2])
            else:
                consensussub= consensustrunksub(entry[1],entry[2])
    converttoyear=True
    if converttoyear:
        humanM1 = convertbranchlengthtoyear('%s/human/M1/trunksubstitutionsummary.csv'% os.getcwd(),'%s/human/M1/trunksubstitutionyearsummary.csv'% os.getcwd())
        humanNP = convertbranchlengthtoyear('%s/human/NP/trunksubstitutionsummary.csv'% os.getcwd(), '%s/human/NP/trunksubstitutionyearsummary.csv'% os.getcwd())
        swineM1 = convertbranchlengthtoyear('%s/swine/M1/trunksubstitutionsummary.csv'% os.getcwd(),'%s/swine/M1/trunksubstitutionyearsummary.csv'% os.getcwd())
        swineNP = convertbranchlengthtoyear('%s/swine/NP/trunksubstitutionsummary.csv'% os.getcwd(),'%s/swine/NP/trunksubstitutionyearsummary.csv'% os.getcwd())

    combfile=True
    if combfile:
        combinedsummaryfile = '%s/plots/S9_trunksubstitutionsummary.csv' % os.getcwd()
        combinefiles(humanM1, humanNP, swineM1,swineNP, combinedsummaryfile)  

        



        


main()