'''
This script performs analysis looking at whether human influenza T-cell epitopes are under selection compared to swine influenza T-cell epitopes.
The statistic we are using to accomplish this we call f. Given toy data as shown below.

Amino acid sites in the protein of interest, the number of unique epitopes at each site (determined using epitopefinder), and
the average number of mutations that occur at that site (determined by averaging values over BEAST trees).

aminoacid_site  number_unique_epitopes  average_human_mutations  average_swine_mutations
1               0                       1                        0
2               3                       2                        1
3               1                       6                        4

f = sum(numberepitopes*mutations)/sum(mutations)
fhuman = (0*1 + 3*2 + 1*3)/(1 + 2 + 6) = 1
fswine = (0*0 + 3*1 + 1*4)/(1 + 1 + 4) = 7/6

We are interested in the difference between human and swine f values. In order to test the significance of this difference,
we create a null distribution by shuffling the number of unique epitopes on each site. This analysis is done for influenza NP and M1 for mutational information derived from the whole tree and for the trunk of the tree.
All of the randomized and actual f values are recorded.

Functions
-----------
*``readfiletodictionary`` : reads files to a dictionary. Used to read in files containing the amino acid sites, number of unique epitopes, and average number of mutatations at each site
*``initiatesummary`` : initiates the output summary file ``randomdistributionsummary.csv``
*``makedirectory`` : checks if directory exists and if it doesn't, the function creates it
*``AvgEpitopeChangesPerSub`` : calculates f
*``ComputeStatAndPValue`` : main function for f analysis. computes f for human and swine 
Writes each value of the null distribution to ``randomdistribution.csv`` and actual f values to ``randomdistributionsummary.csv``

Input files
------------
*``avemutationpersite.csv`` : contains file with sites and average mutations per site (for each protein there is a file for the whole tree and the trunk)
*``combinedepitopesbysite.csv`` : contains file with sites and unique epitopes per site (for each protein)

Output files
--------------
*``randomdistribution.csv`` : file created for whole tree and trunk for each protein. contains null distribution of randomized human f, swine f, and delta f values

'''

import os
import re
import random
import matplotlib

matplotlib.use('pdf')
import pylab
import numpy as np

def readfiletodictionary(filename, makelist = False):
    '''this function takes a ``filename`` where the file is in the following format:
    titleline
    1, 0
    2, 2
    3, 5
    and returns a dictionary with first entry as key and second entry as value
    There is an option to return list containing first column by setting makelist=True when function called
        '''

    fx = open(filename, 'r')
    sites = []
    filedict = {}
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            if makelist:           
                sites.append(int(entry[0]))
            filedict[int(entry[0])] = float(entry[1])
    fx.close()
    if makelist:
        return sites, filedict
    else:
        return filedict

def AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site):
    """Calculates average epitope changes per site.

    *sites* is a list of all sites.

    *epitopes_per_site* is dictionary keyed by each site in *sites* with
    value giving number of epitopes.

    *subs_per_site* is dictionary keyed by each site in *sites* with
    value giving number of epitopes.

    The returned value is the average number of epitopes changed per substitution.

    >>> sites = [1, 2]
    >>> epitopes_per_site = {1:2, 2:0}
    >>> subs_per_site = {1:1, 2:2}
    >>> x = AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site)
    >>> print "%.2f" % x
    0.67
    """
    mtot = weighted_mtot = 0
    for site in sites:
        mtot += subs_per_site[site]
        weighted_mtot += subs_per_site[site] * epitopes_per_site[site]
    if mtot == 0:
        return 0
    return weighted_mtot / float(mtot)


def ComputeStatAndPValue(protein, sites, epitopes_per_site, subs_per_site1, subs_per_site2, nrandom, outputdistribution):
    """This function computes f (through *AvgEpitopeChangesPerSub*) for human and swine and delta f(fhuman - fswine).
    A null distribution of delta f values is created by randomizing the number of epitopes_per_site and recalculating 
    delta f. A one sided  P-value that f for *subs_per_site1* is greater than for *subs_per_site2*. 

    *``protein`` : protein
    *``sites`` : amino acid sites
    *``epitopes_per_site`` : dictionary with sites as key and epitopes per site as the value
    *``subs_per_site1`` : dictionary with sites as key and average substitutions per site for human influenza as the value
    *``subs_per_site2`` : dictionary with sites as key and average substitutions per site for swine influenza as the value
    *``nrandom`` : number of randomizations for the null distribution
    *``outputdistribution`` : name of output file to record all of the f human, f swine, delta f values computed in the randomization

    """
    random_distribution_list = [] 
   
    x1 = AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site1)
    x2 = AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site2)
    actual_delta_x1_x2 = []
   
    random_distribution_list = []

    f = open(outputdistribution, 'w')
    f.write('actual, %s, %s\n' % (x1, x2))
    random.seed(0)
    for irandom in range(nrandom):
        random_epitopes_per_site = list(epitopes_per_site.values())
        random.shuffle(random_epitopes_per_site)
        random_epitopes_per_site_dict = dict(zip(sites, random_epitopes_per_site))

        random_x1 = AvgEpitopeChangesPerSub(sites, random_epitopes_per_site_dict, subs_per_site1)
        random_x2 = AvgEpitopeChangesPerSub(sites, random_epitopes_per_site_dict, subs_per_site2)

        f.write('random, %s, %s\n' % (random_x1, random_x2))  
  
    print " %s %f %f " % (protein, x1, x2)   
    f.close()
    
def initiatesummary(outputsummary):
    '''This function writes the first line of the output summary file
    *outputsummary* summary of f statistic'''

    f = open(outputsummary, 'w')
    f.write('protein, P, humanf, swinef, deltaf\n')
    f.close()

def makedirectory(dirname):
    '''This function checks if a directory exists, and if it doesn't, it creates the directory
    *dirname* directory name 
    '''
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
       

def main():
       
    epitopetypes = ( 'cd8',)
    analysistype = ('trunk','tree')
  
    for atype in analysistype:
        print atype
        masterdistribution = []
        masterdeltaf = []
        masterhumantree = []
        masterswinetree = []
        masterhumantrunk = []
        masterswinetrunk = []
    
        for epitopetype in epitopetypes:
            print epitopetype

            checkdirectory = True
            outputdir = '%s/plots'% (os.getcwd())
            outputsubdir = '%s/plots/%s'%(os.getcwd(), epitopetype)
            #outputsubsubdir = ('%s/plots/%s/%s'%(os.getcwd(),epitopetype,atype))
            if checkdirectory:
                directory = makedirectory(outputdir)
                directory = makedirectory(outputsubdir)
                #directory = makedirectory(outputsubsubdir)

            proteins = ['M1','NP']           
                   
            for protein in proteins:
                humanmutationpersitefile = '%s/human/%s/%savemutationpersite.csv' % (os.getcwd(), protein, atype)
                swinemutationpersitefile = '%s/swine/%s/%savemutationpersite.csv' % (os.getcwd(), protein, atype)
                mutationstodict = True
                if mutationstodict:
                    sites, humanmutpersite = readfiletodictionary(humanmutationpersitefile,makelist=True)
                    swinemutpersite = readfiletodictionary(swinemutationpersitefile)

                epfile = '%s/human/%s/%scombinedepitopesbysite.csv' % (os.getcwd(), protein, epitopetype)
                epitopestodict = True
                if epitopestodict:
                    epitopesbysite = readfiletodictionary(epfile)
                
                computefstat = True
                nulldistributionoutfile = '%s/human/%s/%s%srandomdistribution.csv' % (os.getcwd(), protein, epitopetype,atype)                 
                if computefstat:
                    print 'computing f statistic'
                    numberrandom = int(10000)
                    if atype == 'tree':
                        fstat = ComputeStatAndPValue(protein, sites, epitopesbysite, humanmutpersite, swinemutpersite, numberrandom,nulldistributionoutfile)
                    if atype == 'trunk':
                        fstat= ComputeStatAndPValue(protein, sites, epitopesbysite, humanmutpersite, swinemutpersite, numberrandom,nulldistributionoutfile)





main()
            
