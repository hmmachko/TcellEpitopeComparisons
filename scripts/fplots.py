'''
We are interested if the ratio of fhuman/fswine,ftrunk/ftree and (ftree/ftrunk)human/(ftree/ftrunk)swine is significant for our analysis of selection on 
T-cell epitopes in M1 and NP. In order to test the significance of this difference,we create a null distribution by shuffling the number of unique epitopes 
on each site and calculating a 1-sided p-value. This analysis is donefor fhuman/fswine(trunk), fhuman/fswine(tree), ftrunk/ftree(human), ftrunk/ftree(swine), 
and (ftree/ftrunk)human/(ftree/ftrunk)swine for both M1 and NP. 

We make a scatter plot that displays all  of the f values. For fhuman/fswine, ftrunk/ftree, and (ftree/ftrunk)human/(ftree/ftrunk)swine we make a violin
plot of the null distribution and the actual ratio overlayed as a datapoint. 

Functions
-----------
*``opencsvfile`` : retrieves f values from file
*``ratioandpvalue`` : computes the p value 
*``ratio2`` : computes the ratio of 2 numbers
*``appenddata`` : appends data points to a list
*``appenddata2`` : appends data points to a list
*``makeviolinplot`` : makes violin plots of fhuman/fswine and ftrunk/ftree
*``scatterplot`` : makes plot containing all actual f values
*``summary`` : writes summary of p values for fhuman/fswine and ftrunk/ftree

Input files
-------------
*``cd8treerandomdistribution.csv`` : file containing actual f values and random f values, calculated from whole tree
*``cd8trunkrandomdistribution.csv`` : file containing actual f values and random f values, calculated from trunk

Output files
--------------
*``scatterfvalues.pdf`` : plot of all f values
*``ratiofhumanswine.pdf``: violin plot with null distribution and actual values for fhuman/fswine
*``ratiotrunktree.pdf`` : violin plot with null distribution and actual values for ftree/ftrunk
*``ratiotrunktreehumanswine.pdf`` : violin plot with null distribution and actual values for (ftree/ftrunk)human/(ftree/ftrunk)swine
*``summarypvalues.csv`` : p values for all analyses
'''

import os
import re
import random
import matplotlib
from matplotlib import rc

matplotlib.use('pdf')
import pylab
import numpy as np

def opencsvfile(fname):
    '''This function retrieves data from a ``fname`` csv file and returns the values for the human tree,
    human trunk, swine tree, swine trunk ratio of average epitope change per mutation. 

    *``fname`` : input csv file in a format such that human trunk, human tree, swine trunk, swine tree is in the 
    first entry of a line, the human value is the second entry and the swine value is the third entry.
    '''
    nullhuman = []
    nullswine = []


    fx = open(fname,'r')
    with fx as f:
        #next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            if entry[0] == 'actual':
                actualhuman = float(entry[1])
                actualswine = float(entry[2])
            else:
            
                #protein = entry[0]
                nullhuman.append(float(entry[1]))
                nullswine.append(float(entry[2]))
    fx.close()
    return nullhuman,nullswine,actualhuman,actualswine

def makeviolinplot(datalists,plotfile,masterdeltaf,xlabs,ylab,pos,yax,xtitle):
    '''This function constructs a violin plot using matplotlib with 4 violins paired in groups of 2.
    It also overlays a datapoint on each violin plot. 

    *datalists* list of lists, each list contains the data for a violin plot (the layout works for 
        4 violins that are arranged in groups of 2)
    *plotfile* outputfile of violinplot in pdf format
    *protein* xlabel
    *masterdeltaf* list of 4 points, each placed on a violin plot 
    *xlabs* x labels for pair of plots
    *ylab* y label
    *pos* x axis positions of viiolin plots
    *yax* y axis limits
    *xtitle* title for x axis
    '''

    fig, ax = matplotlib.pyplot.subplots()
    plot = matplotlib.pyplot.violinplot(datalists, pos, points=60, widths=0.4, showmeans=False,
                      showextrema=False, showmedians=True, bw_method="scott")
    plot = matplotlib.pyplot.scatter(pos,masterdeltaf,s=100, color='Red')
    matplotlib.pyplot.xlabel(xtitle,fontsize = 30)
    matplotlib.pyplot.ylim(yax)
    ax.set_xticks(pos)
    ax.set_xticklabels(xlabs,fontsize = 30)
    ax.set_ylabel(ylab,fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 30)
    matplotlib.pyplot.tight_layout()
    if 'testratiotrunktree.pdf' in plotfile:
        matplotlib.pyplot.yticks(np.arange(0, 3.1, 1))

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()  

def scatterplot(data1,data2,data3,data4,data5,data6,data7,data8,plotfile):
    ''' This function creates a scatterplot with some connected lines. For our specific example, we are 
    plotting f values for human and swine trunk and tree for M1 and NP. On the x axis, we plot M1 human, 
    M1 swine, NP human, NP swine. Thus the trunk and tree values are both plotted at the same x-coordinate.
    The values for trunk human and swine (for a protein) have a line connecting them. The same is true for tree 
    human and swine.
    
    *data1* M1actualhumantree
    *data2* M1actualswinetree
    *data3* NPactualhumantree
    *data4* NPactualswinetree
    *data5* M1actualhumantrunk
    *data6* M1actualswinetrunk
    *data7* NPactualhumantrunk
    *data8* NPactualswinetrunk
    *plotfile* pdf name of plot

    '''

    datapointsa = []
    datapointsb = []
    datapointsc = []
    datapointsd = []
    datapointsa.append(data1)
    datapointsa.append(data2)
    datapointsb.append(data3)
    datapointsb.append(data4)
    datapointsc.append(data5)
    datapointsc.append(data6)
    datapointsd.append(data7)
    datapointsd.append(data8)

    x_val = [1,2]
    x_val2 = [3,4]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 

    matplotlib.pyplot.plot(x_val, datapointsa, '-o',label = 'tree',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsb, '-o',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val, datapointsc, '-o',label = 'trunk',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsd, '-o',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.ylabel('$F$',fontsize = 30)
    matplotlib.pyplot.xlabel('M1                    NP',fontsize = 30)
    matplotlib.pyplot.legend(loc=2, fontsize = 30)
    xlabs = ['human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.yticks(fontsize = 30)
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 30)
    matplotlib.pyplot.xlim([0.5, 4.5])
    matplotlib.pyplot.ylim([0, 1.6])
    matplotlib.pyplot.yticks(np.arange(0, 1.7, .4))
    matplotlib.pyplot.tight_layout()
    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()   
    
def ratioandpvalue(data1,data2,actual):
    ratio12 = []
    n=0
    for i,value in enumerate(data1):
        ratio = float(value)/data2[i] 
        ratio12.append(ratio)
        if ratio >= actual:
                n += 1

    p = max(n/float(len(data1)), 1 / float(len(data1)))
    print p

    return ratio12,p

def ratio2(data1,data2):
    '''This function takes the ratio of *data1* to data2* and returns that ratio'''

    ratio12 = float(data1)/data2
    return ratio12

    

def appenddata(data1,data2,data3,data4):
    '''This function appends *data1*, *data2*, *data3*, and *data4* to a list in that order
    and returns the list.'''

    master = []
    master.append(data1)
    master.append(data2)
    master.append(data3)
    master.append(data4)
    return master

def appenddata2(data1,data2):
    '''This function appends *data1* and *data2* to a list in that order and returns the list.'''
    master = []
    master.append(data1)
    master.append(data2)

    return master

def summary(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,outputsummary):
    '''This function writes the first line of the output summary file
    *outputsummary* summary of f statistic'''

    f = open(outputsummary, 'w')
    f.write('summary P values\n')
    f.write(',M1 tree,M1 trunk,NP tree,NP trunk\n')
    f.write('fhuman/fswine,%f,%f,%f,%f\n' %(p1,p2,p3,p4))
    f.write(',M1 human, M1 swine, NP human, NP swine, HA human, HA swine\n')
    f.write('ftrunk/ftree,%f,%f,%f,%f\n' %(p5,p6,p7,p8))
    f.write(',M1,NP\n')
    f.write('fhuman/fswine,%f,%f\n' %(p9,p10))
    f.close()

def main():

    datafiles = (
        '%s/human/M1/cd8treerandomdistribution.csv' % (os.getcwd()),
        '%s/human/M1/cd8trunkrandomdistribution.csv' % (os.getcwd()),
        '%s/human/NP/cd8treerandomdistribution.csv' % (os.getcwd()),
        '%s/human/NP/cd8trunkrandomdistribution.csv' % (os.getcwd()),
    )
        
    readdata = True
    if readdata:
        print 'reading files'
        M1humantree,M1swinetree,M1actualhumantree,M1actualswinetree =  opencsvfile(datafiles[0])
        M1humantrunk,M1swinetrunk,M1actualhumantrunk,M1actualswinetrunk=  opencsvfile(datafiles[1])
        NPhumantree,NPswinetree,NPactualhumantree,NPactualswinetree=  opencsvfile(datafiles[2])
        NPhumantrunk,NPswinetrunk,NPactualhumantrunk,NPactualswinetrunk=  opencsvfile(datafiles[3])
        
    scatter = True
    plotfile = '%s/plots/cd8/scatterfvalues.pdf' % (os.getcwd())
    if scatter:
        print 'making scatterplot'
        splot = scatterplot(M1actualhumantree,M1actualswinetree,NPactualhumantree,NPactualswinetree,M1actualhumantrunk,M1actualswinetrunk,NPactualhumantrunk, NPactualswinetrunk,plotfile)

    calchumanswinefratio = True
    if calchumanswinefratio:

        actualM1tree = ratio2(M1actualhumantree,M1actualswinetree)
        actualM1trunk = ratio2(M1actualhumantrunk,M1actualswinetrunk)
        actualNPtree = ratio2(NPactualhumantree,NPactualswinetree)
        actualNPtrunk = ratio2(NPactualhumantrunk,NPactualswinetrunk)      
        M1tree,pM1tree = ratioandpvalue(M1humantree,M1swinetree,actualM1tree)
        M1trunk,pM1trunk = ratioandpvalue(M1humantrunk,M1swinetrunk,actualM1trunk)
        NPtree,pNPtree = ratioandpvalue(NPhumantree,NPswinetree,actualNPtree)
        NPtrunk,pNPtrunk = ratioandpvalue(NPhumantrunk,NPswinetrunk,actualNPtrunk)
       
    appendfratio = True
    if appendfratio:
        humanswinefratiodata = appenddata(M1trunk,M1tree,NPtrunk,NPtree)
        actualhumswinefractiodata = appenddata(actualM1trunk,actualM1tree,actualNPtrunk,actualNPtree)

    plothumanswinefratio = True
    pos =  [1,2,3.5,4.5]
    xlabs = ['trunk','tree','trunk','tree']
    xtitle = 'M1                     NP'
    ylab = '$F_{human}/F_{swine}$'
    plotfile = '%s/plots/cd8/ratiofhumanswine.pdf' % (os.getcwd())
    yax = [-0.5, 6]
    if plothumanswinefratio:
        print 'plotting fhuman/fswine'
        humanswinefratio = makeviolinplot(humanswinefratiodata,plotfile,actualhumswinefractiodata,xlabs,ylab,pos,yax,xtitle)

    plottrunktotree = True
    plotfile = '%s/plots/cd8/ratiotrunktree.pdf' % (os.getcwd())
    xlabs = ['human','swine','human','swine']
    ylab = '$F_{trunk}/F_{tree}$' #r_{ijk} #'$\frac{F_{trunk}}/{F_{tree}}$'
    yax = [-0.1, 3]
    if plottrunktotree:
        actualM1humantrunktotree = ratio2(M1actualhumantrunk,M1actualhumantree)
        actualM1swinetrunktotree= ratio2(M1actualswinetrunk,M1actualswinetree)
        actualNPhumantrunktotree= ratio2(NPactualhumantrunk,NPactualhumantree)
        actualNPswinetrunktotree= ratio2(NPactualswinetrunk,NPactualswinetree)
        M1humantrunktotree,pM1human = ratioandpvalue(M1humantrunk,M1humantree,actualM1humantrunktotree)
        M1swinetrunktotree,pM1swine = ratioandpvalue(M1swinetrunk,M1swinetree,actualM1swinetrunktotree)
        NPhumantrunktotree,pNPhuman = ratioandpvalue(NPhumantrunk,NPhumantree,actualNPhumantrunktotree)
        NPswinetrunktotree,pNPswine = ratioandpvalue(NPswinetrunk,NPswinetree,actualNPswinetrunktotree)
        trunktreefratiodata = appenddata(M1humantrunktotree,M1swinetrunktotree,NPhumantrunktotree, NPswinetrunktotree)
        actualtrunktreefratiodata = appenddata(actualM1humantrunktotree,actualM1swinetrunktotree,actualNPhumantrunktotree, actualNPswinetrunktotree)
        print 'plotting ftrunk/tree'
        trunktreefratio = makeviolinplot(trunktreefratiodata,plotfile,actualtrunktreefratiodata,xlabs,ylab,pos,yax,xtitle)

    plottrunktotreehumantoswine = True
    plotfile = '%s/plots/cd8/ratiotrunktreehumanswine.pdf' % (os.getcwd())
    xlabs = ['M1','NP']
    xtitle = ''
    ylab = '$(F_{trunk}/F_{tree})_{h}/(F_{trunk}/F_{tree})_{s}$'
    pos =  [1,2]
    yax = [-0.5, 10]
    if plottrunktotreehumantoswine:

        actualM1trunktotreehumantoswine = ratio2(actualM1humantrunktotree,actualM1swinetrunktotree)
        actualNPtrunktotreehumantoswine = ratio2(actualNPhumantrunktotree,actualNPswinetrunktotree)
        M1trunktotreehumantoswine,pM1 = ratioandpvalue(M1humantrunktotree,M1swinetrunktotree,actualM1trunktotreehumantoswine)
        NPtrunktotreehumantoswine,pNP = ratioandpvalue(NPhumantrunktotree,NPswinetrunktotree,actualNPtrunktotreehumantoswine)
        trunktreehumanswinefratiodata = appenddata2(M1trunktotreehumantoswine,NPtrunktotreehumantoswine)
        actualtrunktreehumanswinefratiodata = appenddata2(actualM1trunktotreehumantoswine,actualNPtrunktotreehumantoswine)
        trunktreehumanswinefratio = makeviolinplot(trunktreehumanswinefratiodata,plotfile,actualtrunktreehumanswinefratiodata,xlabs,ylab,pos,yax,xtitle)
        count = 0
        countM1 = 0
        for value in NPtrunktotreehumantoswine:
            if value > actualNPtrunktotreehumantoswine:
                count +=1
        for value in M1trunktotreehumantoswine:
            if value > actualM1trunktotreehumantoswine:
                countM1 +=1

        print count, countM1

    writesummary = True
    if writesummary:
        summarypvalues = '%s/plots/cd8/summarypvalues.csv' % (os.getcwd())
        pvalsummary = summary(pM1tree,pM1trunk,pNPtree,pNPtrunk,pM1human,pM1swine,pNPhuman,pNPswine,pM1,pNP,summarypvalues)














    #nulldistributionoutfile = '%s/human/%s/%s%srandomdistribution.csv' % (os.getcwd(), protein, epitopetype,atype)


main()