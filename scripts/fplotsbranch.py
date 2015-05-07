#fplotting

import os
import re
import random
import matplotlib
from matplotlib import rc

matplotlib.use('pdf')
import pylab
import numpy as np

#read in output
#actual x1,x2
#random x1,x2

def opencsvfile(fname):
    '''This function retrieves data from a ``fname`` csv file and returns the values for the human tree,
    human trunk, swine tree, swine trunk ratio of average mutation per epitope codon vs nonepitope codon. 

    *``fname`` : input csv file in a format such that human trunk, human tree, swine trunk, swine tree is in the 
    first entry of a line and the desired value is in the second entry
    *``order`` : list containing the desired order for returning data values from file 
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
    '''This function constructs a violin plot using matplotlib

    *datalists* lists of datapoints
    *plotfile* outputfile of violinplot in pdf format
    *protein* xlabel
    *masterdeltaf* list of points, 1 placed on each violin plot 

    '''
   # #pos =  [1,2,3]
    #fig = matplotlib.pyplot.figure()
    #ax = fig.add_subplot(111)
    
    #bp = ax.boxplot(data,0,'')

    fig, ax = matplotlib.pyplot.subplots()
    plot = matplotlib.pyplot.violinplot(datalists, pos, points=60, widths=0.4, showmeans=False,
                      showextrema=False, showmedians=True, bw_method="scott")
    #plot = matplotlib.pyplot.boxplot(datalists, notch=0, sym='', positions=pos,vert=1, whis=1.5)
    plot = matplotlib.pyplot.scatter(pos,masterdeltaf,s=100, color='Red')
    matplotlib.pyplot.xlabel(xtitle,fontsize = 22)
    #rc('text', usetex=True)
    matplotlib.pyplot.ylim(yax)
    ax.set_xticks(pos)
    ax.set_xticklabels(xlabs,fontsize = 22)
    ax.set_ylabel(ylab,fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 18)
    matplotlib.pyplot.tight_layout()
    #ax.gca.set_ylabel(r'$\boldsymbol{\delta}average epitope change per mutation$')

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()  

def scatterplot(data1,data2,data3,data4,data5,data6,data7,data8,plotfile):
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
    datapointsc = []
    datapointsd = []
    datapointse =[]
    datapointsf = []
    #datapointsa.append(data1)
    #datapointsa.append(data2)
    #datapointsa.append(data3)
    #datapointsa.append(data4)
   # datapointsb.append(data5)
    #datapointsb.append(data6)
   # datapointsb.append(data7)
   # datapointsb.append(data8)

    datapointsa.append(data1)
    datapointsa.append(data2)
    datapointsb.append(data3)
    datapointsb.append(data4)
    datapointsc.append(data5)
    datapointsc.append(data6)
    datapointsd.append(data7)
    datapointsd.append(data8)
    #datapointse.append(data9)
    #datapointse.append(data10)
    #datapointsf.append(data11)
    #datapointsf.append(data12)

    #x_val = [1,2,4,5]
    x_val = [1,2]
    x_val2 = [3,4]
    x_val3 = [5,6]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 

    matplotlib.pyplot.plot(x_val, datapointsa, '-o',label = 'tree',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsb, '-o',color = [210/256.0,94/256.0,202/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val, datapointsc, '-o',label = 'trunk',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, datapointsd, '-o',color = [94/256.0,210/256.0,114/256.0],markersize=15,linewidth=3)
   # matplotlib.pyplot.plot(x_val3, datapointse, '-o',color = [95/256.0,158/256.0,209/256.0])
   # matplotlib.pyplot.plot(x_val3, datapointsf, '-o',color = [255/256.0,128/256.0,14/256.0])

    #matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'tree')
    #matplotlib.pyplot.scatter(x_val, datapointsb,color = [255/256.0,128/256.0,14/256.0], label = 'trunk')
    matplotlib.pyplot.ylabel('$F$ (ave epitope change per mutation)',fontsize = 22)
    matplotlib.pyplot.xlabel('M1                              NP',fontsize = 22)
    matplotlib.pyplot.legend(loc=2, fontsize = 20)
    xlabs = ['human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.yticks(fontsize = 18)
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 22)
    matplotlib.pyplot.xlim([0.5, 4.5])
    matplotlib.pyplot.ylim([0, 1.6])
    matplotlib.pyplot.yticks(np.arange(0, 1.7, .4))
    matplotlib.pyplot.tight_layout()


    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()   

def deltahumanswine(data1,data2):
    ratio12 = []
    for i,value in enumerate(data1):
        ratio = float(value)/data2[i] 
        ratio12.append(ratio)

    return ratio12

def ratio(data1,data2,actual):
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

    ratio12 = float(data1)/data2
    return ratio12
def actualdeltahumanswine(data1,data2):
    ratio12 = float(data1)/data2
    return ratio12

def appenddata(data1,data2,data3,data4):
    master = []
    master.append(data1)
    master.append(data2)
    master.append(data3)
    master.append(data4)
    #master.append(data5)
    #master.append(data6)
    return master

def appenddata2(data1,data2):
    master = []
    master.append(data1)
    master.append(data2)
    #master.append(data3)

    return master

def summary(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,outputsummary):
    '''This function writes the first line of the output summary file
    *outputsummary* summary of f statistic'''

    f = open(outputsummary, 'w')
    f.write('summary P values\n')
    f.write(',M1 branch,M1 trunk,NP branch,NP trunk,HA branch,HA trunk\n')
    f.write('fhuman/fswine,%f,%f,%f,%f\n' %(p1,p2,p3,p4))
    f.write(',M1 human, M1 swine, NP human, NP swine, HA human, HA swine\n')
    f.write('ftrunk/ftree,%f,%f,%f,%f\n' %(p5,p6,p7,p8))
    f.write(',M1,NP,HA\n')
    f.write('fhuman/fswine,%f,%f\n' %(p9,p10))
    f.close()


def main():

    datafiles = (
        '%s/human/M1/cd8branchrandomdistributionbranch.csv' % (os.getcwd()),
        '%s/human/M1/cd8trunkrandomdistributionbranch.csv' % (os.getcwd()),
        '%s/human/NP/cd8branchrandomdistributionbranch.csv' % (os.getcwd()),
        '%s/human/NP/cd8trunkrandomdistributionbranch.csv' % (os.getcwd()),
     
    )
        
    readdata = True
    if readdata:
        print 'reading files'
        M1humantree,M1swinetree,M1actualhumantree,M1actualswinetree =  opencsvfile(datafiles[0])
        M1humantrunk,M1swinetrunk,M1actualhumantrunk,M1actualswinetrunk=  opencsvfile(datafiles[1])
        NPhumantree,NPswinetree,NPactualhumantree,NPactualswinetree=  opencsvfile(datafiles[2])
        NPhumantrunk,NPswinetrunk,NPactualhumantrunk,NPactualswinetrunk=  opencsvfile(datafiles[3])
       

    scatter = True
    plotfile = '%s/plots/cd8/testscatterbranch.pdf' % (os.getcwd())
    if scatter:
        print 'making scatterplot'
        splot = scatterplot(M1actualhumantree,M1actualswinetree,NPactualhumantree,NPactualswinetree,M1actualhumantrunk,M1actualswinetrunk,NPactualhumantrunk, NPactualswinetrunk,plotfile)


    calchumanswinefratio = True
    if calchumanswinefratio:

        actualM1tree = actualdeltahumanswine(M1actualhumantree,M1actualswinetree)
        actualM1trunk = actualdeltahumanswine(M1actualhumantrunk,M1actualswinetrunk)
        actualNPtree = actualdeltahumanswine(NPactualhumantree,NPactualswinetree)
        actualNPtrunk = actualdeltahumanswine(NPactualhumantrunk,NPactualswinetrunk)
        
        M1tree,pM1tree = ratio(M1humantree,M1swinetree,actualM1tree)
        M1trunk,pM1trunk = ratio(M1humantrunk,M1swinetrunk,actualM1trunk)
        NPtree,pNPtree = ratio(NPhumantree,NPswinetree,actualNPtree)
        NPtrunk,pNPtrunk = ratio(NPhumantrunk,NPswinetrunk,actualNPtrunk)
       

    appendfratio = True
    if appendfratio:
        humanswinefratiodata = appenddata(M1trunk,M1tree,NPtrunk,NPtree)
        actualhumswinefractiodata = appenddata(actualM1trunk,actualM1tree,actualNPtrunk,actualNPtree)

    plothumanswinefratio = True
    pos =  [1,2,4,5]
    xlabs = ['trunk','branch','trunk','branch']
    xtitle = 'M1                              NP'
    ylab = '$F_{human}/F_{swine}$'
    plotfile = '%s/plots/cd8/testratiofhumanswinebranch.pdf' % (os.getcwd())
    yax = [-0.5, 8]
    if plothumanswinefratio:
        print 'plotting fhuman/fswine'
        humanswinefratio = makeviolinplot(humanswinefratiodata,plotfile,actualhumswinefractiodata,xlabs,ylab,pos,yax,xtitle)

    plottrunktotree = True
    plotfile = '%s/plots/cd8/testratiotrunkbranch.pdf' % (os.getcwd())
    xlabs = ['human','swine','human','swine']
    ylab = '$F_{trunk}/F_{branch}$' #r_{ijk} #'$\frac{F_{trunk}}/{F_{tree}}$'
    yax = [-0.1, 5]
    if plottrunktotree:
        actualM1humantrunktotree = ratio2(M1actualhumantrunk,M1actualhumantree)
        actualM1swinetrunktotree= ratio2(M1actualswinetrunk,M1actualswinetree)
        actualNPhumantrunktotree= ratio2(NPactualhumantrunk,NPactualhumantree)
        actualNPswinetrunktotree= ratio2(NPactualswinetrunk,NPactualswinetree)
        
        M1humantrunktotree,pM1human = ratio(M1humantrunk,M1humantree,actualM1humantrunktotree)
        M1swinetrunktotree,pM1swine = ratio(M1swinetrunk,M1swinetree,actualM1swinetrunktotree)
        NPhumantrunktotree,pNPhuman = ratio(NPhumantrunk,NPhumantree,actualNPhumantrunktotree)
        NPswinetrunktotree,pNPswine = ratio(NPswinetrunk,NPswinetree,actualNPswinetrunktotree)
        
        trunktreefratiodata = appenddata(M1humantrunktotree,M1swinetrunktotree,NPhumantrunktotree, NPswinetrunktotree)
        actualtrunktreefratiodata = appenddata(actualM1humantrunktotree,actualM1swinetrunktotree,actualNPhumantrunktotree, actualNPswinetrunktotree)
        #print len(trunktreefratiodata[0]),len(trunktreefratiodata[1]),len(trunktreefratiodata[2]),len(trunktreefratiodata[3])
       # print trunktreefratiodata[0]
        #print trunktreefratiodata[1]
       # print trunktreefratiodata[2]
       # print trunktreefratiodata[3]
       # print len(trunktreefratiodata)
       # print actualtrunktreefratiodata
        print 'plotting ftrunk/tree'
        trunktreefratio = makeviolinplot(trunktreefratiodata,plotfile,actualtrunktreefratiodata,xlabs,ylab,pos,yax,xtitle)

    plottrunktotreehumantoswine = True
    plotfile = '%s/plots/cd8/testratiotrunkbranchhumanswine.pdf' % (os.getcwd())
    xlabs = ['M1','NP']
    xtitle = ''
    ylab = '$(F_{trunk}/F_{branch})_{h}/(F_{trunk}/F_{branch})_{s}$'
    pos =  [1,2]
    yax = [-0.5, 10]
    if plottrunktotreehumantoswine:

        actualM1trunktotreehumantoswine = ratio2(actualM1humantrunktotree,actualM1swinetrunktotree)
        actualNPtrunktotreehumantoswine = ratio2(actualNPhumantrunktotree,actualNPswinetrunktotree)
        
        M1trunktotreehumantoswine,pM1 = ratio(M1humantrunktotree,M1swinetrunktotree,actualM1trunktotreehumantoswine)
        NPtrunktotreehumantoswine,pNP = ratio(NPhumantrunktotree,NPswinetrunktotree,actualNPtrunktotreehumantoswine)
        
        trunktreehumanswinefratiodata = appenddata2(M1trunktotreehumantoswine,NPtrunktotreehumantoswine)
        actualtrunktreehumanswinefratiodata = appenddata2(actualM1trunktotreehumantoswine,actualNPtrunktotreehumantoswine)
        trunktreehumanswinefratio = makeviolinplot(trunktreehumanswinefratiodata,plotfile,actualtrunktreehumanswinefratiodata,xlabs,ylab,pos,yax,xtitle)
        #trunktreehumanswinefratio = makeviolinplot(NPtrunktotreehumantoswine,plotfile,actualNPtrunktotreehumantoswine,xlabs,ylab,pos)
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
        summarypvalues = '%s/plots/cd8/summarypvaluesbranch.csv' % (os.getcwd())
        pvalsummary = summary(pM1tree,pM1trunk,pNPtree,pNPtrunk,pM1human,pM1swine,pNPhuman,pNPswine,pM1,pNP,summarypvalues)














    #nulldistributionoutfile = '%s/human/%s/%s%srandomdistribution.csv' % (os.getcwd(), protein, epitopetype,atype)


main()