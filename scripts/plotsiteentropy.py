
import os
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
    count = 0
    epitopecount = 0
    fx = open(filename, 'r')
    sites = []
    filedict = {}
    with fx as f:
        next(f)
        for lines in fx: 
            count+=1    
            entry = lines.strip().split(",")
            #print entry
            if makelist:           
                sites.append(int(entry[0]))
            if int(entry[1]) != 0:
                epitopecount+=1
            filedict[int(entry[0])] = int(entry[1])
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
            print entry


            filedict[int(entry[0])] = float(entry[2])
    fx.close()
    return filedict

def sortbynumberepitopes(epitopedict,siteentropydict,summaryfile):
    from collections import defaultdict
    averageentropydict = {}
    entropybyepitopenumber = {0:[[],[]],1:[[],[]],2:[[],[]],3:[[],[]],4:[[],[]],5:[[],[]],6:[[],[]],7:[[],[]],8:[[],[]],9:[[],[]],}
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

    for epnumber in entropybyepitopenumber:
        if epnumber !=0:
            for entry in entropybyepitopenumber[epnumber][1]:
                totent_ep += entry
                count_ep +=1
    average_ep = float(totent_ep)/count_ep
    print 'average for nonepitopes %f' % averageentropydict[0]
    print 'average for epitopes %f' % average_ep
    f = open(summaryfile,'w')
    for epnumber in sorted(averageentropydict):
        f.write('%s,%s\n' % (epnumber,averageentropydict[epnumber]))
    f.write('average epitope, %s' % average_ep)


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
        print len(siteentropydict[entry][0]), len(siteentropydict[entry][1])
       # print siteentropydict[entry][0]
       # print siteentropydict[entry][1]
        matplotlib.pyplot.scatter(siteentropydict[entry][0],siteentropydict[entry][1])
        #print ok
    count=0
    for entry in aveentropydict:
        matplotlib.pyplot.scatter(count, aveentropydict[count],color='red')
        count+=1
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
    matplotlib.pyplot.ylim(-1,6)
    matplotlib.pyplot.yticks(np.arange(0, 6, 1))
    matplotlib.pyplot.tight_layout()
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
main()



