'''
The purpose of this script is to create a plot of the ratio of substitution rate in epitope regions to
non-epitope regions. This analysis is completed for influenza M1 and NP for both data 
for the entire tree and the trunk of the tree.

Functions:
------------
*``opencsvfile`` : this function retrieves the previously calculated ratio of the average substitution rate in epitopes to nonepitopes
*``scatterplot`` : this function creates a summary graph of the average substitution rate in epitopes vs non-epitopes

Input files:
-------------
*``avemutationratio_epitope_nonepitope.csv`` : contains previously calculated data from BEAST trees. 

Output files:
--------------
*``epitope_to_nonepitope_mutation.pdf`` : summary graph of the average substitution rate in epitopes vs non-epitopes

Written by Heather Machkovech March 2015
'''

import os
import re
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np

def opencsvfile(fname):
    '''This function retrieves data from a *fname* csv file and returns the values for the human tree,
    human trunk, swine tree, swine trunk ratio of average mutation per epitope codon vs nonepitope codon. 

    *``fname`` : input csv file in a format such that human trunk, human tree, swine trunk, swine tree is in the 
    first entry of a line and the desired value is in the second entry
    *``order`` : list containing the desired order for returning data values from file 
    '''
 
    proteins = []
    ratedict = {}
    #ratetrunkdict ={}
    datalines = ['M1 human trunk','M1 swine trunk','NP human trunk','NP swine trunk','M1 human tree','M1 swine tree','NP human tree','NP swine tree']
    fx = open(fname,'r')
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            if entry[0] in datalines:
                ratedict[entry[0]] = (float(entry[1]))
    fx.close()

    M1trunk = []
    M1tree = []
    NPtrunk = []
    NPtree = []
    M1trunk.append(ratedict['M1 human trunk'])
    M1trunk.append(ratedict['M1 swine trunk'])
    M1tree.append(ratedict['M1 human tree'])
    M1tree.append(ratedict['M1 swine tree'])
    NPtrunk.append(ratedict['NP human trunk'])
    NPtrunk.append(ratedict['NP swine trunk'])
    NPtree.append(ratedict['NP human tree'])
    NPtree.append(ratedict['NP swine tree'])

 
    return  M1trunk,M1tree,NPtrunk,NPtree

def scatterplot(datapointsa,datapointsb,datapointsd,datapointse,plotfile):
    '''This function creates a scatterplot of epitope/nonepitope substitution rate
    and saves the pdf to *plotfile*. The input data is *datapointsa*, *datapointsb*, *datapointsd*, and *datapointse*, 
    lists (M1tree, NPtree, M1trunk, NPtrunk)that contain 2 values (human and swine). The plot has a line 
    between the human and swine values.
    '''
    

    x_val = [1,2]
    x_val2 = [3,4]
    
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 
    matplotlib.pyplot.plot(x_val, datapointsa, '-o',label = 'tree',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsb, '-o',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val, datapointsd, '-o',label = 'trunk',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointse, '-o',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.ylabel('epitope/nonepitope rate',fontsize = 30)
    matplotlib.pyplot.xlabel('M1                    NP',fontsize = 30)
    matplotlib.pyplot.legend(loc=4, fontsize = 30)
    xlabs = ['human','swine','human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 30)
    matplotlib.pyplot.axhline(y=1,linewidth=2,color='Black')
    matplotlib.pyplot.xlim([0.5, 4.5])
    matplotlib.pyplot.ylim([-.1, 1.1])
    matplotlib.pyplot.yticks(fontsize = 30)
    matplotlib.pyplot.yticks(np.arange(0, 1.1, .25))
    matplotlib.pyplot.tight_layout()

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()   

def main():
     
    infile1 = '%s/plots/cd8/avemutationratio_epitope_nonepitope.csv' % (os.getcwd())
    outfile1 = '%s/plots/cd8/epitope_to_nonepitope_mutation_ratio.pdf' %(os.getcwd())
      
    yax = 'epitope to nonepitope mutation ratio'
    readcsv = True
    if readcsv:        
        M1trunk,M1tree,NPtrunk,NPtree = opencsvfile(infile1)
          
    bar = True
    if bar:
        graph1 = scatterplot(M1tree,NPtree,M1trunk,NPtrunk,outfile1)


main()