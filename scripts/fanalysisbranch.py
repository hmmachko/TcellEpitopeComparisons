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
we create a null distribution by shuffling the number of unique epitopes on each site and calculating a 
1-sided p-value. This analysis is done for influenza NP, PA, and M1 for mutational information derived from the whole tree and for the trunk of the tree.

For each protein, a histogram plot is constructed that contains the null distribution with an overlay of the actual delta f and the p-value.
2 summary graphs are created. One is a paired bar graph of human and swine f values for each protein. The second is a  violin plot is 
constructed that contains a violin for each protein (representing the null distrubution of delta f) and an additional scatter point of the 
actual deltaf is overlayed. 

Functions
-----------
*``readfiletodictionary`` : reads files to a dictionary. Used to read in files containing the amino acid sites, number of unique epitopes, and average number of mutatations at each site
*``initiatesummary`` : initiates the output summary file ``randomdistributionsummary.csv``
*``makedirectory`` : checks if directory exists and if it doesn't, the function creates them
*``AvgEpitopeChangesPerSub`` : calculates f
*``ComputeStatAndPValue`` : main function for f analysis. computes f for human and swine, delta f, null distribution, and 1-sided p-value. 
Writes each value of the null distribution to ``randomdistribution.csv`` and summary of 
*``histogramplot`` : creates a histogram plot for each protein containing the null delta f distribution, a line representing the actual delta f, and the p-value
*``makeviolinplot`` : creates a violin plot with a violin for each protein representing the delta f null distribution, with an additional scatter point of the actual delta f overlayed
*``pairedbargraph`` : function that generates paired bar graph for human and swine f values

Input files
------------
*``avemutationpersite.csv`` : contains file with sites and average mutations per site (for each protein there is a file for the whole tree and the trunk)
*``combinedepitopesbysite.csv`` : contains file with sites and unique epitopes per site (for each protein)

Output files
--------------
*``randomdistributionsummary.csv`` : summary file created for whole tree and trunk. contains: f human, f swine, delta f, and p-value
*``randomdistribution.csv`` : file created for whole tree and trunk for each protein. contains null distribution of randomized human f, swine f, and delta f values
*``distribution.pdf`` : histogram plot of null distribution for whole tree and trunk for each protein. 
*``violindistribution.pdf`` : a violin plot containing all proteins, 1 for whole tree and 1 for trunk
*``fhuman_vs_fswine.pdf`` : bar graph of human and swine f values, 1 for whole tree and 1 for trunk


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


def ComputeStatAndPValue(protein, sites, epitopes_per_site, subs_per_site1, subs_per_site2, nrandom, outputdistribution, outputsummary,masterdistribution,masterdeltaf,masterhuman,masterswine):
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
    *``outputsummary`` : name of summary output file containing f human, f swine, delta f, and p-value
    *``masterdistribution`` : list of all values of null distribution for a protein, used for graphing
    *``masterdeltaf`` : list of all the delta f values for all proteins, used for graphing
    *``masterhuman`` : list of all the human f values for all proteins, used for graphing
    *``masterswine`` : list of all the swine f values for all proteins, used for graphing

    """
    random_distribution_list = []
    norm_random_distribution_list = []

    n_ge = 0
   
    x1 = AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site1)
    x2 = AvgEpitopeChangesPerSub(sites, epitopes_per_site, subs_per_site2)
    actual_delta_x1_x2 = []
   
    delta_x1_x2 = x1 - x2 
    
    masterhuman.append(x1)
    masterswine.append(x2)
   
    actual_delta_x1_x2.append(delta_x1_x2)
    random_distribution_list = []


    f = open(outputdistribution, 'w')
    fx = open(outputsummary, 'a')
    fx.write('%s,' % protein)
    f.write('actual, %s, %s, %s\n' % (x1, x2, delta_x1_x2))
    random.seed(0)
    for irandom in range(nrandom):
        random_epitopes_per_site = list(epitopes_per_site.values())
        random.shuffle(random_epitopes_per_site)
        random_epitopes_per_site_dict = dict(zip(sites, random_epitopes_per_site))

        random_x1 = AvgEpitopeChangesPerSub(sites, random_epitopes_per_site_dict, subs_per_site1)
        random_x2 = AvgEpitopeChangesPerSub(sites, random_epitopes_per_site_dict, subs_per_site2)

        random_delta_x1_x2 = random_x1 - random_x2
        random_distribution_list.append(random_delta_x1_x2)
        if (random_delta_x1_x2 >= delta_x1_x2):
                n_ge += 1
        f.write('random, %s, %s, %s \n' % (random_x1, random_x2, random_delta_x1_x2))  
        if irandom%1000== 0:
            print random_x1,random_x2   

    p = max(n_ge / float(nrandom), 1 / float(nrandom))
  
    print " %s %f %f %f %f" % (protein, x1, x2, delta_x1_x2,p)
    fx.write('%s,%s,%s,%s\n' % (p,x1, x2, delta_x1_x2))
   
    f.close()
    fx.close()
    masterdistribution.append(random_distribution_list)
    masterdeltaf.append(actual_delta_x1_x2)

    return (actual_delta_x1_x2, random_distribution_list, p,masterhuman,masterswine)

def histogramplot(datalists, dataline, plotfile,protein,pval): 
    '''This function creates a histogram plot. 
    *datalists* list of values for histogram 
    *dataline* value for vertical line to display on graph (in this case the actual f value)
    *plotfile* ouput file for plot
    *protein* used for title
    *pval* p-value displayed on graph
    '''
   
    prot = re.sub('_', ' ',protein)
    pylab.hist(datalists,color = [95/256.0,158/256.0,209/256.0]) #color = 'DodgerBlue'
    pylab.figtext(0.8,0.85,"p=%s"%pval)
    pylab.axvline(dataline,color = 'Black')
    pylab.xlabel('average epitope change per mutation',fontsize = 15)
    pylab.ylabel('number of randomizations',fontsize = 15)
    #pylab.title('%s' % (prot))
    matplotlib.pyplot.yticks(fontsize = 11)
    matplotlib.pyplot.xticks(fontsize = 11)


    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

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

def makeviolinplot(datalists,plotfile,protein,masterdeltaf,xlabs,ylab):
    '''This function constructs a violin plot using matplotlib

    *datalists* lists of datapoints
    *plotfile* outputfile of violinplot in pdf format
    *protein* xlabel
    *masterdeltaf* list of points, 1 placed on each violin plot 

    '''
    pos =  [1,2,3]
    fig, ax = matplotlib.pyplot.subplots()
    plot = matplotlib.pyplot.violinplot(datalists, pos, points=60, widths=0.4, showmeans=False,
                      showextrema=False, showmedians=False, bw_method="scott")
    plot = matplotlib.pyplot.scatter(pos,masterdeltaf,s=60, color='Red')
    ax.set_xticks(pos)
    ax.set_xticklabels(xlabs)
    ax.set_ylabel(ylab)
    #ax.gca.set_ylabel(r'$\boldsymbol{\delta}average epitope change per mutation$')

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()        

def pairedbargraph(data, data2,data3,data4, plotfile, yaxis,xlabs):

    '''This function creates a paired bar graph.

    *data* is a list containing the first data points of the pair 
    *data2* is a list containing the second data points of the pair
    *plotfile* is the pdf filename for the plot
    *yaxis* title for yaxis
    *xlabs* labels for x ticks
    '''
    x = xrange(len(data))
    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.2
    print data
    print data2
    print data3
    print data4

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, color = [0,107/256.0,164/256.0],label = 'human tree') #color = "DarkBlue"
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width, color = [200.0/256,82/256.0,0],label = 'swine tree')#color = "DarkRed"

    cplot = matplotlib.pyplot.bar(index + 2*bar_width,data3, bar_width, color = [95/256.0,158/256.0,209/256.0],label = 'human trunk')#color = "DarkCyan"
    dplot = matplotlib.pyplot.bar(index + 3*bar_width,data4, bar_width,color = [255/256.0,128/256.0,14/256.0],label = 'swine trunk')#color = "Chocolate"


    matplotlib.pyplot.ylabel(yaxis,fontsize = 15)
    matplotlib.pyplot.yticks(fontsize = 11)
    matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)
    matplotlib.pyplot.legend(loc=2, fontsize = 13)
   # matplotlib.pyplot.axhline(linewidth=4)        # inc. width of x-axis 
    #matplotlib.pyplot.axvline(linewidth=4)

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def scatterplot(data1,data2,data3,data4,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk
    #xticks - human/swine
    #M1: 
       # human: 1,data1(tree)
               # 1,data2(trunk)
       # swine: 2,data1(tree)
               # 2,data2(trunk)
    #NP:
       # human: 4,data1(tree)
               # 4,data2(trunk)
        #swine: 5,data1(tree)
               # 5,data2(trunk)

    datapointsa = []
    datapointsb = []
    datapointsa.append(data1[0])
    datapointsa.append(data3[0])
    datapointsa.append(data1[1])
    datapointsa.append(data3[1])
    datapointsb.append(data2[0])
    datapointsb.append(data4[0])
    datapointsb.append(data2[1])
    datapointsb.append(data4[1])
    x_val = [1,2,4,5]
    matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'tree')
    matplotlib.pyplot.scatter(x_val, datapointsb,color = [255/256.0,128/256.0,14/256.0], label = 'trunk')
    matplotlib.pyplot.ylabel('f',fontsize = 15)
    matplotlib.pyplot.legend(loc=2, fontsize = 13)
    xlabs = ['human','swine','human','swine']
    xtick = [1,2,4,5]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 15)


    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)


    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def plotftreevstrunk(data1,data2,data3,data4,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk
    datapointsa = []
    datapointsb = []
    datapointsa.append(float(data2[0])/data1[0])
    datapointsa.append(float(data2[1])/data1[1])
    datapointsb.append(float(data4[0])/data3[0])
    datapointsb.append(float(data4[1])/data3[1])

    x = xrange(len(datapointsa))
    n_groups = len(datapointsa)
    index = np.arange(n_groups)
    bar_width = 0.35

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,datapointsa, bar_width, color = [0,107/256.0,164/256.0],label = 'human') #color = "DarkBlue"
    bplot = matplotlib.pyplot.bar(index+ bar_width,datapointsb, bar_width, color = [200.0/256,82/256.0,0],label = 'swine')#color = "DarkRed"
    matplotlib.pyplot.legend(loc=2, fontsize = 13)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

    return datapointsa,datapointsb



def plotfhumantreetrunkswinetreetrunk(dataa,datab,plotfile):
    values = []
    values.append(float(dataa[0])/(datab[0]))
    values.append(float(dataa[1])/(datab[1]))
    x = xrange(len(values))
    n_groups = len(values)
    index = np.arange(n_groups)
    bar_width = 0.35

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,values, bar_width, color = [0,107/256.0,164/256.0])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()



def opencsvfileb(fname):
    '''This function retrieves data from a ``fname`` csv file and returns the values for the human tree,
    human trunk, swine tree, swine trunk ratio of average mutation per epitope codon vs nonepitope codon. 

    *``fname`` : input csv file in a format such that human trunk, human tree, swine trunk, swine tree is in the 
    first entry of a line and the desired value is in the second entry
    *``order`` : list containing the desired order for returning data values from file 
    '''
    human = []
    swine = []

    fx = open(fname,'r')
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            
                #protein = entry[0]
            human.append(float(entry[2]))
            swine.append(float(entry[3]))
    fx.close()
    return human,swine

def main():
    
    epitopetypes = ('cd8',)
    #epitopetypes = ( 'cd8',)
    analysistype = ('trunk','branch')

  
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
            outputsubsubdir = ('%s/plots/%s/%s'%(os.getcwd(),epitopetype,atype))
            if checkdirectory:
                directory = makedirectory(outputdir)
                directory = makedirectory(outputsubdir)
                directory = makedirectory(outputsubsubdir)

            if epitopetype == 'antibody':
                proteins = ('HA_H3',)
            else:
                proteins = ['M1','NP']   

            initializesummary = True          
            outputsummary = '%s/plots/%s/%s/randomdistributionsummarybranch.csv' % (os.getcwd(), epitopetype, atype)
            if initializesummary:
                summary = initiatesummary(outputsummary)
          
            for protein in proteins:
                humanmutationpersitefile = '%s/human/%s/%savemutationpersite.csv' % (os.getcwd(), protein, atype)
                swinemutationpersitefile = '%s/swine/%s/%savemutationpersite.csv' % (os.getcwd(), protein, atype)
                mutationstodict = True
                if mutationstodict:
                    sites, humanmutpersite = readfiletodictionary(humanmutationpersitefile,makelist=True)
                    swinemutpersite = readfiletodictionary(swinemutationpersitefile)

                if epitopetype == 'antibody':
                    epfile = '%s/human/%s/antibodyepitopesbysite.csv' % (os.getcwd(), protein)
                else:
                    epfile = '%s/human/%s/%scombinedepitopesbysite.csv' % (os.getcwd(), protein, epitopetype)
                epitopestodict = True
                if epitopestodict:
                    epitopesbysite = readfiletodictionary(epfile)
                
                computefstat = True
                nulldistributionoutfile = '%s/human/%s/%s%srandomdistributionbranch.csv' % (os.getcwd(), protein, epitopetype,atype)                 
                if computefstat:
                    print 'computing f statistic'
                    numberrandom = int(10000)
                    if atype == 'branch':
                        actuallist, randomlist,pvalue, masterhumantree,masterswinetree = ComputeStatAndPValue(protein, sites, epitopesbysite, humanmutpersite, swinemutpersite, numberrandom,nulldistributionoutfile,outputsummary,masterdistribution,masterdeltaf,masterhumantree,masterswinetree)
                    if atype == 'trunk':
                        actuallist, randomlist,pvalue, masterhumantrunk,masterswinetrunk= ComputeStatAndPValue(protein, sites, epitopesbysite, humanmutpersite, swinemutpersite, numberrandom,nulldistributionoutfile,outputsummary,masterdistribution,masterdeltaf,masterhumantrunk,masterswinetrunk)

                plothistogram = False
                plotfile = '%s/plots/%s/%s/%sdistributionbranch.pdf' % (os.getcwd(), epitopetype, atype,protein)

                if plothistogram:
                    datalists = [randomlist]
                    makeplot = histogramplot(datalists, actuallist,plotfile,protein,pvalue)

        violinplot = False
        vplotfile = '%s/plots/%s/%s/violindistribution.pdf' % (os.getcwd(), epitopetype, atype)
        ylab = 'epitope change per mut human - epitope change per mut swine'
        if violinplot:
            vplot = makeviolinplot(masterdistribution, vplotfile,proteins,masterdeltaf,proteins,ylab)

    graphfhumanvsfswine = False
    pairedbarfile = '%s/plots/%s/fhuman_vs_fswine.pdf' % (os.getcwd(), epitopetype)

    readcsv = False
    if readcsv:
        humantrunk,swinetrunk = opencsvfileb('%s/plots/cd8/trunk/randomdistributionsummary.csv'% (os.getcwd()))
        humantree, swinetree = opencsvfileb('%s/plots/cd8/tree/randomdistributionsummary.csv'% (os.getcwd()))


    yaxis = 'average epitope change per mutation'
    if graphfhumanvsfswine:
        fhumanvsfswine = pairedbargraph(humantree, swinetree, humantrunk,swinetrunk,pairedbarfile, yaxis,proteins)

    scatter = False
    plotfile = '%s/plots/%s/testscatter.pdf' % (os.getcwd(), epitopetype)
    if scatter:
        splot = scatterplot(humantree, swinetree, humantrunk,swinetrunk,plotfile)

    ratiotrunktree = False
    plotfile = '%s/plots/%s/testtrunktree.pdf' % (os.getcwd(), epitopetype)
    if ratiotrunktree:
        humantrunktree,swinetrunktree = plotftreevstrunk(humantree, swinetree, humantrunk,swinetrunk,plotfile)

    ratiotrunktreehumanswine = False
    plotfile = '%s/plots/%s/testtrunktreehumanswine.pdf' % (os.getcwd(), epitopetype)
    if ratiotrunktreehumanswine:
        graph = plotfhumantreetrunkswinetreetrunk(humantrunktree,swinetrunktree,plotfile)




main()
            
