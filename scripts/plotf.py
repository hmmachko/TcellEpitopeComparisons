'''
This script generates a paired bar graph of human and swine f values for influenza proteins M1, PA, and NP. 
This script creates 2 graphs; one for data including the whole tree and the other for the trunk of the tree. 

Functions
------------
*``opencsvfile`` : function that reads file containing f values for human and swine and returns a separate list of f human values and f swine values
*``pairedbargraph`` : function that generates paired bar graph for human and swine f values

Input files
-------------
*``randomdistributionsummary.csv`` : file containing f values, a separate file for both the whole tree and trunk of the tree is used

Output files
--------------
*``fvalues.pdf`` : paired bar graph, a separate file for both the whole tree and trunk of the tree

Written by Heather Machkovech Feb 2015
'''

import os
import re
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np
import math

def opencsvfile(filename,order):
    '''This function returns 2 lists (f values for humans and f values for swine) obtained from *filename* 
    *``filename`` : file which has a title line and the f value for humans as the 3rd entry and for swine as the 4th entry
    *``order`` : list of proteins for which f values are retrieved
    '''
    fhuman = {}
    fswine= {}

    humanf = []
    swinef= []

    fx = open(filename, 'r')
    
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            fhuman[(entry[0])] = (entry[2])
            fswine[(entry[0])] = (entry[3])
    fx.close()

    for protein in order:
        humanf.append(fhuman[protein])
        swinef.append(fswine[protein])
            
    return humanf,swinef

def pairedbargraph(data, data2,data3,data4,labs, plotfile, yaxis,title):
    '''This function creates a paired bar graph.
    *data* is a list containing the first data points of the pair 
    *data2* is a list containing the second data points of the pair
    *plotfile* is the pdf filename for the plot
    *yaxis* title for yaxis
    *title* graph title
    '''
    
    x = xrange(len(data))

    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.2
    
    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, color = "DodgerBlue", label = 'human tree')
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width,color = "DarkRed", label = 'swine tree')

    cplot = matplotlib.pyplot.bar(index + 2*bar_width,data3, bar_width, color = "DarkCyan",label = 'human trunk')
    dplot = matplotlib.pyplot.bar(index + 3*bar_width,data4, bar_width, color = "Chocolate",label = 'swine trunk')


    matplotlib.pyplot.ylabel(yaxis,fontsize = 15, fontweight='bold')
    matplotlib.pyplot.yticks(fontsize = 11, fontweight='bold')
    matplotlib.pyplot.xticks(index + 2*bar_width, labs, fontsize = 15, fontweight='bold')
    matplotlib.pyplot.legend(fontsize = 13)
    matplotlib.pyplot.axhline(linewidth=4)        # inc. width of x-axis 
    matplotlib.pyplot..axvline(linewidth=4)

    matplotlib.pyplot.title(title)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()




def main():
    analysistype = ('tree', 'trunk')
    epitopetypes = ('cd8',)
    order = ('PA', 'M1','NP')    
    yax = 'f'

    
    for epitope in epitopetypes:
        outfile = ('%s/plots/%s/fvalues.pdf' %(os.getcwd(),epitope))
        filenametree = ('%s/plots/%s/tree/randomdistributionsummary.csv' %(os.getcwd(),epitope))
        filenametrunk = ('%s/plots/%s/trunk/randomdistributionsummary.csv' %(os.getcwd(),epitope))


        readfile = True
        if readfile:
            if epitope =='cd8':
                fvalhtree, fvalstree = opencsvfile(filenametree,order)
                fvalhtrunk, fvalstrunk = opencsvfile(filenametrunk,order)


            barplot = True
            if barplot:
                if epitope =='cd8':          
                    graph1 = pairedbargraph(fvalhtree, fvalstree, fvalhtrunk, fvalstrunk,order, outfile, yax,epitope)
               

main()











