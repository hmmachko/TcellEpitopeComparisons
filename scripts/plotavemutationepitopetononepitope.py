'''
The purpose of this script is to create a bar plot the ratio of the average number of mutations 
per epitope codon to nonepitope codon to compare human and swine influenza values. 

This analysis is completed for influenza PA, M1, NP for both data for the entire tree and the trunk of the tree.

Functions:
------------
*``opencsvfile`` : this function retrieves the previously calculated values of the average number of mutations in epitope codons and nonepitope codons
*``bargraph`` : this function creates a paired bar graph (human vs swine) of the ratio of the average number of mutations 
per epitope codon to nonepitope codon

Input files:
-------------
*``avemutationratio_epitope_nonepitope.csv`` : contains previously calculated data from BEAST trees. 
There are separate files for tree and trunk that contain human and swine values for each protein. 

Output files:
--------------
*``epitope_to_nonepitope_mutation.pdf`` : paired bar graph (human vs swine) of the ratio of the average number of mutations 
per epitope codon to nonepitope codon. Separate graphs for tree and trunk data are generated.

Written by Heather Machkovech March 2015
'''

import os
import re
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np

def opencsvfile(fname):
    '''This function retrieves data from a ``fname`` csv file and returns the values for the human tree,
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

def openantibodyfile(fname):
    proteins = []
    ratedict = {}
    datalines2 = ['HA_H3 human trunk','HA_H3 swine trunk','HA_H3 human tree','HA_H3 swine tree']
    fx = open(fname,'r')
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            if entry[0] in datalines2:
                ratedict[entry[0]] = (float(entry[1]))
    fx.close()

    HAtrunk = []
    HAtree = []
    HAtrunk.append(ratedict['HA_H3 human trunk'])
    HAtrunk.append(ratedict['HA_H3 swine trunk'])
    HAtree.append(ratedict['HA_H3 human tree'])
    HAtree.append(ratedict['HA_H3 swine tree'])

    return  HAtrunk,HAtree

def bargraph(data, data2,data3,data4,labs, plotfile, yaxis):
    '''This function creates a paired bar graph using matplotlib.

    *``data`` : list containing first set of bars
    *``data2`` : list containing second set of bars
    *``labs`` : list of labels for each bar pair
    *``plotfile`` : file name for pdf of plot
    *``yaxis`` : label for y axis
    
    '''
   
    x = xrange(len(data))

    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.2

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, color=[210/256.0,94/256.0,202/256.0],label = 'human tree') #color = "DarkBlue"
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width,color = [94.0/256,82/256.0,0],label = 'swine tree') #color = "DarkRed"
    cplot = matplotlib.pyplot.bar(index + 2*bar_width,data3, bar_width,color = [95/256.0,158/256.0,209/256.0],label = 'human trunk') #color = "DarkCyan"
    dplot = matplotlib.pyplot.bar(index + 3*bar_width,data4, bar_width,color = [255/256.0,128/256.0,14/256.0],label = 'swine trunk') #color = "Chocolate"

    #matplotlib.pyplot.xlabel('protein')
    matplotlib.pyplot.ylabel(yaxis,fontsize = 20)
    #matplotlib.pyplot.title(epitopetype)
    #plot.tick_params(axis='both', which='major', labelsize=10)
    matplotlib.pyplot.yticks(fontsize = 18)

    matplotlib.pyplot.xticks(index + 2*bar_width, labs, fontsize = 20)
    matplotlib.pyplot.legend(fontsize = 18)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])
    matplotlib.pyplot.axhline(y=1,linewidth=2, color = 'Black')


    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def dataappend(masterhumantree, masterhumantrunk, masterswinetree, masterswinetrunk,humantreevalues,swinetreevalues,humantrunkvalues,swinetrunkvalues):
    masterhumantree.append(humantreevalues)    
    masterhumantrunk.append(humantrunkvalues)
    masterswinetrunk.append(swinetrunkvalues)
    masterswinetree.append(swinetreevalues)

    return masterhumantree,masterhumantrunk,masterswinetree, masterswinetrunk

def scatterplot(datapointsa,datapointsb,datapointsc,datapointsd,datapointse,datapointsf,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk


    #x_val = [1,2,4,5]
    x_val = [1,2]
    x_val2 = [3,4]
    x_val3 = [5,6]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 
    matplotlib.pyplot.plot(x_val, datapointsa, '-o',label = 'tree',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsb, '-o',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    #matplotlib.pyplot.plot(x_val3, datapointsc, '-o',color = [95/256.0,158/256.0,209/256.0])
    matplotlib.pyplot.plot(x_val, datapointsd, '-o',label = 'trunk',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointse, '-o',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
   # matplotlib.pyplot.plot(x_val3, datapointsf, '-o',color = [255/256.0,128/256.0,14/256.0])
    #matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'tree')
    #matplotlib.pyplot.scatter(x_val, datapointsb,color = [255/256.0,128/256.0,14/256.0], label = 'trunk')
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
    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()   

def main():
     
    infile1 = '%s/plots/cd8/avemutationratio_epitope_nonepitope.csv' % (os.getcwd())
    infile2 = '%s/plots/antibody/avemutationratio_epitope_nonepitope.csv' % (os.getcwd())
            #infile2 = ('%s/plots/%s/%s/hsrandomdistributionsummary.csv' % (os.getcwd(), epitope, atype))

    outfile1 = '%s/plots/cd8/epitope_to_nonepitope_mutation_ratio.pdf' %(os.getcwd())
        #outfile2 = ('%s/plots/%s/trunkepitope_to_nonepitope_mutation.pdf' %(os.getcwd(),epitope))


    yax = 'epitope to nonepitope mutation ratio'
    readcsv = True
    if readcsv:
            
        M1trunk,M1tree,NPtrunk,NPtree = opencsvfile(infile1)
            
        HAtrunk,HAtree = openantibodyfile(infile2)

    appenddata = False
    if appenddata:
         masterhumantree,masterhumantrunk,masterswinetree, masterswinetrunk = dataappend(masterhumantree, masterhumantrunk, masterswinetree, masterswinetrunk,humantreevalues,swinetreevalues,humantrunkvalues,swinetrunkvalues)


           
    bar = True
    if bar:
        #if epitope == 'cd8':
        graph1 = scatterplot(M1tree,NPtree,HAtree,M1trunk,NPtrunk,HAtrunk,outfile1)
            #graph2 = bargraph(humantrunkvalues, swinetrunkvalues, order1, outfile2, yax)


main()