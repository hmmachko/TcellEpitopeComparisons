'''This script creates a maximum clade credibility tree from a BEAST trees
 file using the BEAST treeannotator. It assumes a 10% burnin 
Input files 
------------
*``prot_aligned.trees`` : BEAST file containing thinned trees 

Output files
--------------
*``maxcladecredibility.tre`` : maximum clade credibility tree
'''
import os
import sys

def main():
    
    proteinlineage_directories = (
    #'%s/human/M1/prot_aligned.trees' % os.getcwd(),  
   # '%s/human/NP/prot_aligned.trees' % os.getcwd(), 
    #'%s/swine/M1/prot_aligned.trees' % os.getcwd(),  
    #'%s/swine/NP/prot_aligned.trees' % os.getcwd(),  
    '%s/human/HA_H3/prot_aligned.trees' % os.getcwd(),  
    '%s/swine/HA_H3/prot_aligned.trees' % os.getcwd(), 
    
    )

    for proteinlineage in proteinlineage_directories:
        makemaxtree = True
        base = proteinlineage[:-18]
        home = os.path.expanduser("~")
        treepath = '%s/BEASTv1.8.1/bin/treeannotator' % home
        if makemaxtree:
            makemaxcladecredtree = os.system('%s -burnin 3000000 %s %smaxcladecredibility.tre' % (treepath, proteinlineage, base))

main()


