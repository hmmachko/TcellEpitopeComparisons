'''This script is to build a quick RAxML tree of combined human and swine lineage to remove any potential outliers.
This analysis was done for influenza  M1 and NP. 

Functions
----------
*``raxmltree`` : makes a RAxML tree of aligned DNA sequences

Input files
------------
*``cds_aligned.fasta`` : aligned DNA sequences (separate human and swine files for each protein)

Output files
--------------
*``combined_cds_aligned.fasta`` : combined, aligned human and swine influenza sequences used to build tree
*``prelimTree`` : RAxML tree

Written by Heather Machkovech Feb 2015

'''

import subprocess
import os
import glob
def raxmltree(sequencesets):
    """ This function generates a quick RAxML tree to visualize anomalous sequences that do not fit the molecular clock.
    sequencesets = fasta file containing sequences for RAxML tree

    """
    use_existing_output = False
    basename = sequencesets[:-26]
    raxml_dir = '%sRAxML_output/combined/' % basename
    treefile = 'prelimTree'
    raxml_tree = "%sRAxML_bestTree.%s" % (raxml_dir, treefile)
    if not os.path.isdir(raxml_dir):
        os.mkdir(raxml_dir)
    print "Using RAxML to build a quick phylogenetic tree for visual analysis of potentially anomalous sequences..." 
    print "\nFor sequence set %s" % sequencesets
    if use_existing_output and os.path.isfile(raxml_tree):
        print "The existing tree of %s will be used, and no new tree will be created." % raxml_tree
    else:
        oldfiles = glob.glob("%s/*%s*" % (raxml_dir, treefile))
        for fname in oldfiles:
            os.remove(fname)
        os.system('%s -w %s -n %s -p 1 -m GTRCAT -s %s --no-bfgs' % ('raxmlHPC-SSE3', raxml_dir, treefile, sequencesets))
    if not os.path.isfile(raxml_tree):
        raise ValueError("Failed to generated expected output tree file %s" % raxml_tree)
    print "RAxML tree complete"


def main():

    humanseqfiles = (
        '%s/human/M1/cds_aligned.fasta' % (os.getcwd()),
        '%s/human/NP/cds_aligned.fasta' % (os.getcwd()),       
    )

    swineseqfiles = (
        '%s/swine/M1/cds_aligned.fasta' % (os.getcwd()),
        '%s/swine/NP/cds_aligned.fasta' % (os.getcwd()),     
    )

    combinedfiles = list(zip(humanseqfiles,swineseqfiles))
    concat = True
    combinedtree = True
    if concat:
        for combinedfile in combinedfiles:
            base = combinedfile[0][:-17]
            combined = os.system('cat %s %s > %scombined_cds_aligned.fasta' % (combinedfile[0], combinedfile[1],base))

            if combinedtree:
                tree = raxmltree('%scombined_cds_aligned.fasta'% base)



main()

        