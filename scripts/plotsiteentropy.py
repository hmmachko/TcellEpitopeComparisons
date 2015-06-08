
import os
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats import spearmanr
import scipy.stats as ss

def readfiletodictionary(filename, makelist = False):
    '''this function takes a ``filename`` where the file is in the following format:
    titleline
    1, 0
    2, 2
    3, 5
    and returns a dictionary with first entry as key and second entry as value
    There is an option to return list containing first column by setting makelist=True when function called
        '''
    count = 0
    epitopecount = 0
    fx = open(filename, 'r')
    sites = []
    filedict = {}
    with fx as f:
        next(f)
        for lines in fx: 
            #count+=1    
            entry = lines.strip().split(",")
            #print entry
            if makelist:           
                sites.append(int(entry[0]))
            if int(entry[1]) != 0:
                epitopecount+=1
            if int(entry[0]) != 1:
                filedict[int(entry[0])] = int(entry[1])
                count+=1 
    fx.close()
    if makelist:
        return sites, filedict
    else:
        return filedict, count,epitopecount

def readfiletodictionary2(filename):
    '''this function takes a ``filename`` where the file is in the following format:
    titleline
    1, 0
    2, 2
    3, 5
    and returns a dictionary with first entry as key and second entry as value
    There is an option to return list containing first column by setting makelist=True when function called
        '''
  
    fx = open(filename, 'r')

    filedict = {}
    with fx as f:
        next(f)
        for lines in fx: 
              
            entry = lines.strip().split()
            #print entry
            if int(entry[0]) != 1:

                filedict[int(entry[0])] = float(entry[2])
    fx.close()
    return filedict

def sortbynumberepitopes(epitopedict,siteentropydict,summaryfile):
    from collections import defaultdict
    averageentropydict = {}
    entropybyepitopenumber = {0:[[],[],[]],1:[[],[],[]],2:[[],[],[]],3:[[],[],[]],4:[[],[],[]],5:[[],[],[]],6:[[],[],[]],7:[[],[],[]],8:[[],[],[]],9:[[],[],[]]}
    for site in epitopedict:
        entropybyepitopenumber[epitopedict[site]][0].append(epitopedict[site])
        entropybyepitopenumber[epitopedict[site]][1].append(siteentropydict[site])
    for epnumber in entropybyepitopenumber:
        totalentropy = 0
        for entry in entropybyepitopenumber[epnumber][1]:
            totalentropy += entry
        averageent = float(totalentropy)/len(entropybyepitopenumber[epnumber][1])
        averageentropydict[epnumber] = averageent
    totent_ep=0
    count_ep =0
    entropyforrank = []
    #ranksforentropy = []
    for epnumber in sorted(entropybyepitopenumber):
        for sentropy in entropybyepitopenumber[epnumber][1]:
            entropyforrank.append(sentropy)
    ranksforentropy = ss.rankdata(entropyforrank)
    print ranksforentropy

    for epnumber in sorted(entropybyepitopenumber):
        for sentropy in entropybyepitopenumber[epnumber][1]:
            indexforrank = entropyforrank.index(sentropy)
            entropybyepitopenumber[epnumber][2].append(ranksforentropy[indexforrank])


    for epnumber in entropybyepitopenumber:
        if epnumber !=0:
            for entry in entropybyepitopenumber[epnumber][1]:
                totent_ep += entry
                count_ep +=1
    average_ep = float(totent_ep)/count_ep
   # print 'average for nonepitopes %f' % averageentropydict[0]
   # print 'average for epitopes %f' % average_ep
    f = open(summaryfile,'w')
    for epnumber in sorted(averageentropydict):
        f.write('%s,%s\n' % (epnumber,averageentropydict[epnumber]))
    f.write('average epitope, %s\n' % average_ep)
    f.close()


    #for epnumber in entropybyepitopenumber:
        #print entropybyepitopenumber[epnumber][0]
        #print entropybyepitopenumber[epnumber][1]
    return entropybyepitopenumber, averageentropydict

def scatterplot(siteentropydict,plotfile,aveentropydict):
    ''' This function creates a scatterplot with some connected lines. For our specific example, we are 
    plotting proportion dN/dS values greater than 1for human and swine epitope and nonepitope regions for M1 and NP. 
    On the x axis, we plot M1 human, M1 swine, NP human, NP swine. Thus the epitope and nonepitope values are both plotted at the same x-coordinate.
    The values for epitope human and swine (for a protein) have a line connecting them. The same is true for nonepitope 
    human and swine.
    
    *``datapointsa`` : M1epitopeproportion
    *``datapointsb`` : NPepitopeproportion
    *``datapointsc`` : M1nonepitopeproportion
    *``datapointsd`` : NPnonepitopeproportion
    *plotfile* pdf name of plot

    '''
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk
    x=0
    xo=10
    y=0
    yo=6
    x,xo,y,yo = matplotlib.pyplot.axis()
    plot_margin=.25

   
    #x_val = [1,2]
    #x_val2 = [3,4]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 
    for entry in siteentropydict:
        #print len(siteentropydict[entry][0]), len(siteentropydict[entry][1])
       # print siteentropydict[entry][0]
       # print siteentropydict[entry][1]
        matplotlib.pyplot.scatter(siteentropydict[entry][0],siteentropydict[entry][2],alpha = 0.2,edgecolor='')
        #print ok
    #count=0
    #for entry in aveentropydict:
        #matplotlib.pyplot.scatter(count, aveentropydict[count],color='red')
       # count+=1
    #matplotlib.pyplot.plot(x_val2, b, '-o',color = [95/256.0,158/256.0,209/256.0],markersize=15,linewidth=3)
    #matplotlib.pyplot.plot(x_val, c, '-o',label = 'nonepitope',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)
    #matplotlib.pyplot.plot(x_val2, e, '-o',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)

    #print 'values for plotting'
    #print datapointsa,datapointsb,datapointsc,datapointsd
    matplotlib.pyplot.ylabel('site entropy',fontsize = 30)
    matplotlib.pyplot.xlabel('number of epitopes',fontsize = 30)
    #if 'FEL' in plotfile:
       # matplotlib.pyplot.legend(loc=4, fontsize = 30)
    #else:
       # matplotlib.pyplot.legend(fontsize = 30)
    #xlabs = ['human','swine','human','swine','human','swine']
    xtick = [0,1,2,3,4,5,6,7,8,9]
    matplotlib.pyplot.xticks(xtick, fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 30)   
    matplotlib.pyplot.xlim([-0.5, 9.5])
    matplotlib.pyplot.ylim(0,500)
    matplotlib.pyplot.yticks(np.arange(0, 501, 100))
    matplotlib.pyplot.tight_layout()
    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close() 

def listsforcorr(dicta,dictb):
    lista = []
    listb = []
    for site in sorted(dicta):
        lista.append(dicta[site])
        listb.append(dictb[site])

    return lista,listb


def pearsoncorr(lista, listb,summaryfile):
    pcorrcoef = pearsonr(lista,listb)
    print 'pearson %s %s' % (pcorrcoef[0],pcorrcoef[1])
    
    f=open(summaryfile,'a')
    f.write('pearsoncorrelation,p-value\n')
    f.write('%s,%s\n'% (pcorrcoef[0],pcorrcoef[1]))
    f.close()
def spearmancorr(lista, listb,summaryfile):
    scorrcoef = spearmanr(lista,listb)
    print 'spearman %s %s' % (scorrcoef[0],scorrcoef[1])
    f=open(summaryfile,'a')
    f.write('spearmancorrelation,p-value\n')
    f.write('%s,%s\n'% (scorrcoef[0],scorrcoef[1]))
    f.close()

def makeviolinplot(datadict,plotfile,xlabs,ylab,pos,yax,xtitle):
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
    datalists = []
    for epitopenumber in sorted(datadict):
        if epitopenumber !=9:
            datalists.append(datadict[epitopenumber][1])
    fig, ax = matplotlib.pyplot.subplots()
    plot = matplotlib.pyplot.violinplot(datalists, pos, points=60, widths=0.4, showmeans=False,
                      showextrema=False, showmedians=True, bw_method="scott")
    plot = matplotlib.pyplot.scatter(9,datadict[9][1],s=100, color='Red')
    matplotlib.pyplot.xlabel(xtitle,fontsize = 30)
    matplotlib.pyplot.ylim(yax)
    matplotlib.pyplot.xlim(-0.5,9.5)
    ax.set_xticks(xlabs)
    ax.set_xticklabels(xlabs,fontsize = 30)
    ax.set_ylabel(ylab,fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 30)
    matplotlib.pyplot.tight_layout()
    
    matplotlib.pyplot.yticks(np.arange(0, 4.5, 1))

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()
def main():
    epitopefiles ='%s/human/NP/cd8combinedepitopesbysite.csv' % os.getcwd()

    siteentropyfile = '%s/inputfiles/meansiteentropy.csv' % os.getcwd()
    entropyplot = '%s/plots/cd8/entropybyepitopenumber.pdf'% os.getcwd()
    summaryentropy = '%s/plots/cd8/entropysummary.csv'% os.getcwd()
    getepitopes = True
    if getepitopes:
        NPepitopes,NPproteinlength,NPepitopesites = readfiletodictionary(epitopefiles, makelist = False)


    getentropies = True
    if getentropies:
        siteentropy = readfiletodictionary2(siteentropyfile)

    sortbyepitopenumber=True
    if sortbyepitopenumber:
        sortedentropies,averageentropy = sortbynumberepitopes(NPepitopes,siteentropy,summaryentropy)

    plotentropyep=True
    if plotentropyep:
        entropyplotbyep = scatterplot(sortedentropies,entropyplot,averageentropy)

    pcorr = True
    if pcorr:
        epitopelist, entropylist = listsforcorr(NPepitopes, siteentropy)
        pcc = pearsoncorr(epitopelist, entropylist,summaryentropy)
        scc = spearmancorr(epitopelist, entropylist,summaryentropy)
        #xlabs= [0,1,2,3,4,5,6,7,8,9]
        #ylab = 'site entropy'
        #pos = [0,1,2,3,4,5,6,7,8]
        #yax = [0,4.5]
        #xtitle = 'number of epitopes'
        #plotentropyviolin = makeviolinplot(sortedentropies,entropyplot,xlabs,ylab,pos,yax,xtitle)
main()



