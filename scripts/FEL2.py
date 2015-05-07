
import os
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np
import sys

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

def separateepitopeandnonepitopevalues(filename,epitopedict,outfile,proteinlength,epitopecount,epitopeproportion,nonepitopeproportion):
    
    fx = open(filename, 'r')
    dnds_epitope = []
    dnds_nonepitope = []
    norm_dnds_epitope = []
    norm_dnds_nonepitope = []
    print 'reading file'
    for lines in fx:
        lines = lines.replace('\r','\n')
        #print lines
    f = open(outfile,'w')
    f.write(lines)
    f.close()

    print 'checking adjusted file'
    epzero = 0
    nonepzero = 0
    epundefined = 0
    nonepundefined = 0
    epgreaterthan1 = 0
    nonepgreaterthan1 = 0
    f = open(outfile,'r')
    with f as fx:
        next(fx)
        for lines in fx:
            entry = lines.strip().split(',')
            #print entry
            #print entry[3]
            if epitopedict[int(entry[0])] != 0:
                if entry[3] != 'Undefined':
                    if entry[3] == 'Infinite':
                        epgreaterthan1 +=1
                    elif float(entry[3]) > 1:
                        epgreaterthan1 +=1
                    else:
                        continue

               # if entry[3] == 'Undefined':
                   #epundefined+=1  
                    #norm_dnds_epitope.append(float(1))
             
               # elif entry[3] =='0':
                   # norm_dnds_epitope.append(float(entry[3]))
   
                   # epzero +=1
               # elif (entry[3] == 'Infinite'):
                   # norm_dnds_epitope.append(float(sys.maxint))               

                #else:
                   # norm_dnds_epitope.append(float(entry[3]))
            else:
                if entry[3] != 'Undefined':

                    if entry[3] == 'Infinite':
                        nonepgreaterthan1 +=1
                    elif float(entry[3]) > 1:
                        nonepgreaterthan1 +=1
                    else:
                        continue
                #if entry[3] == 'Undefined':
                   # nonepundefined+=1  
                    #norm_dnds_nonepitope.append(float(1))
             
               # elif entry[3] =='0':
                    #nonepzero +=1
                   # norm_dnds_nonepitope.append(float(entry[3]))
               # elif entry[3] == 'Infinite':

                 #   norm_dnds_nonepitope.append(float(sys.maxint))
               # else:
                   # norm_dnds_nonepitope.append(float(entry[3]))

   # print 'proportion 0 for epitope %f' % (float(epzero)/epitopecount)
   # print 'proportion 0 for nonepitope %f' % (float(nonepzero)/(proteinlength -epitopecount))
    #print 'proportion undefined for epitope %f' % (float(epundefined)/epitopecount)
    #print 'proportion undefined for nonepitope %f' % (float(nonepundefined)/(proteinlength -epitopecount))

    proportiongreaterthan1epitope = float(epgreaterthan1)/epitopecount
    proportiongreaterthan1nonepitope = float(nonepgreaterthan1)/(proteinlength -epitopecount)
    epitopeproportion.append(proportiongreaterthan1epitope)
    nonepitopeproportion.append(proportiongreaterthan1nonepitope)
  
    totalep = 0
    totalnonep = 0
    for entry in norm_dnds_epitope:
        totalep+=entry
    for entry in norm_dnds_nonepitope:
        totalnonep+= entry


   # print 'epitope'
   # print norm_dnds_epitope 
   # print totalep
   # print float(totalep)/len(norm_dnds_epitope)
   # print 'nonepitope'
   # print norm_dnds_nonepitope 
   # print totalnonep
   # print float(totalnonep)/len(norm_dnds_nonepitope)

    return epitopeproportion, nonepitopeproportion


   # return norm_dnds_epitope, norm_dnds_nonepitope
    #return data

def boxplot(data,labels,outfile):
    print 'plotting'
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    
    bp = ax.boxplot(data,0,'') 

    ax.set_xticklabels(labels)
    pylab.show()
    pylab.savefig(outfile)
    pylab.clf()
    pylab.close()  

def barplot(data,labels,outfile):
    print 'plotting'
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    
    bp = ax.bar(data,0,'') 

    ax.set_xticklabels(labels)
    pylab.show()
    pylab.savefig(outfile)
    pylab.clf()
    pylab.close()  

def bargraph(data, data2,labs, plotfile, yaxis):
   
    x = xrange(len(data))

    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.35

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, color = "DodgerBlue", label = 'epitope')
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width,color = "FireBrick", label = 'nonepitope')
    #matplotlib.pyplot.xlabel('protein')
    matplotlib.pyplot.ylabel(yaxis)
    #matplotlib.pyplot.title(epitopetype)
    matplotlib.pyplot.xticks(index + bar_width, labs)
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def dataappend(masterdata,human_n_dnds_epitope,human_n_dnds_nonepitope,swine_n_dnds_epitope,swine_n_dnds_nonepitope):
    masterdata.append(human_n_dnds_epitope)    
    masterdata.append(human_n_dnds_nonepitope)
    masterdata.append(swine_n_dnds_epitope)
    masterdata.append(swine_n_dnds_nonepitope)

    return masterdata
def makedirectory(dirname):
    '''This function checks if a directory exists, and if it doesn't, it creates the directory
    *dirname* directory name 
    '''
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

def bargraph4(data, data2,data3,data4,labs, plotfile, yaxis):
   
    x = xrange(len(data))

    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.2

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, color = [0,107/256.0,164/256.0],label = 'human epitope')#color = "DodgerBlue"
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width, color = [200.0/256,82/256.0,0], label = 'human nonepitope')#color = "FireBrick"
    cplot = matplotlib.pyplot.bar(index + 2*bar_width,data3, bar_width,color = [95/256.0,158/256.0,209/256.0],label = 'swine epitope') #color = "DarkCyan"
    dplot = matplotlib.pyplot.bar(index + 3*bar_width,data4, bar_width,color = [255/256.0,128/256.0,14/256.0],label = 'swine nonepitope')
    #matplotlib.pyplot.xlabel('protein')
    matplotlib.pyplot.ylabel(yaxis,fontsize = 15)
    #matplotlib.pyplot.title(epitopetype)
    matplotlib.pyplot.xticks(index + 2*bar_width, labs,fontsize = 15)
    matplotlib.pyplot.yticks(fontsize = 11)
    matplotlib.pyplot.legend(loc = 2,fontsize = 13)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def appenddata2(dataa,datab,mastera,masterb,masterc,masterd):
   

    mastera.append(dataa[0])
    masterb.append(datab[0])
    masterc.append(dataa[1])
    masterd.append(datab[1])
    

    return mastera, masterb,masterc,masterd

def scatterplot(datapointsa,datapointsb,datapointsc,datapointsd,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk
    print 'plotting'

    print plotfile
    #x_val = [1,2,4,5]
    x_val = [1,2]
    x_val2 = [3,4]
    matplotlib.pyplot.plot(x_val, datapointsa, '-o',label = 'epitope',color = [95/256.0,158/256.0,209/256.0])
    matplotlib.pyplot.plot(x_val2, datapointsb, '-o',label = 'epitope',color = [95/256.0,158/256.0,209/256.0])
    matplotlib.pyplot.plot(x_val, datapointsc, '-o',label = 'nonepitope',color = [255/256.0,128/256.0,14/256.0])
    matplotlib.pyplot.plot(x_val2, datapointsd, '-o',label = 'nonepitope',color = [255/256.0,128/256.0,14/256.0])

    #matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'tree')
    #matplotlib.pyplot.scatter(x_val, datapointsb,color = [255/256.0,128/256.0,14/256.0], label = 'trunk')
    matplotlib.pyplot.ylabel('fraction of sites with dnds > 1',fontsize = 15)
    matplotlib.pyplot.xlabel('M1                                           NP',fontsize = 15)
    matplotlib.pyplot.legend(loc=1, fontsize = 13)
    xlabs = ['human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 15)
    matplotlib.pyplot.xlim([0.5, 4.5])

    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close() 



def main():

    humanFELfiles = (
        #'%s/dnds/FEL/human_PA_FEL_report.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/human_M1_FEL_report.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/human_NP_FEL_report.csv' % os.getcwd(),
        )
    swineFELfiles = (
        #'%s/dnds/FEL/swine_PA_FEL_report.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/swine_M1_FEL_report.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/swine_NP_FEL_report.csv' % os.getcwd(),
        )
    epitopefiles = (
       # '%s/human/PA/cd8combinedepitopesbysite.csv' % os.getcwd(),
        '%s/human/M1/cd8combinedepitopesbysite.csv' % os.getcwd(),
        '%s/human/NP/cd8combinedepitopesbysite.csv' % os.getcwd(),
        )
    boxplotoutfiles = (
       # '%s/plots/dnds/PA_FELboxplot3.pdf' % os.getcwd(),
        '%s/plots/dnds/M1_FELboxplot3.pdf' % os.getcwd(),
        '%s/plots/dnds/NP_FELboxplot3.pdf' % os.getcwd(),
        )
    adjustedhumanFELfiles = (
       # '%s/dnds/FEL/human_PA_FEL_report_fix.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/human_M1_FEL_report_fix.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/human_NP_FEL_report_fix.csv' % os.getcwd(),
        )
    adjustedswineFELfiles = (
       # '%s/dnds/FEL/swine_PA_FEL_report_fix.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/swine_M1_FEL_report_fix.csv' % os.getcwd(),
        '%s/dnds/FEL/usertree/swine_NP_FEL_report_fix.csv' % os.getcwd(),
        )
    proportionsiteswithhighdndsplotoutfile = (
       # '%s/plots/dnds/PA_highdnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertree/M1_highdnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertree/NP_highdnds.pdf' % os.getcwd(),
        )
    masterproteinproportion_humanepitope = []
    masterproteinproportion_humannonepitope = []
    masterproteinproportion_swineepitope = []
    masterproteinproportion_swinenonepitope = []   
        

    analysis = list(zip(humanFELfiles,swineFELfiles,epitopefiles,boxplotoutfiles,adjustedhumanFELfiles,adjustedswineFELfiles,proportionsiteswithhighdndsplotoutfile))
    #masterdata = []
   # labels = ['human epitope','human nonepitope', 'swine epitope', 'swine nonepitope']
    labs = ['human', 'swine']
    print 'beginning'

    outputdir = ('%s/plots/dnds'%(os.getcwd()))
    checkdirectory = True
    if checkdirectory:
        directory = makedirectory(outputdir)

    for protein in analysis:
        masterdata = []
        epitopeproportion = []
        nonepitopeproportion = []
        readepitopes = True
        if readepitopes:
            epitopes,proteinlength,epitopesites = readfiletodictionary(protein[2], makelist = False)
        getdndsvalues = True
        if getdndsvalues:
            dnds = readfiletodictionary(protein[0])
            epitopeproportion, nonepitopeproportion = separateepitopeandnonepitopevalues(protein[0],epitopes,protein[4],proteinlength,epitopesites,epitopeproportion,nonepitopeproportion)
            epitopeproportion, nonepitopeproportion = separateepitopeandnonepitopevalues(protein[1],epitopes,protein[5],proteinlength,epitopesites,epitopeproportion,nonepitopeproportion)
        appenddata = False
        if appenddata:
            masterdata = dataappend(masterdata, human_n_dnds_epitope,human_n_dnds_nonepitope,swine_n_dnds_epitope,swine_n_dnds_nonepitope)
        makebarplot = False
        yaxis = ('fraction of sites with dnds > 1')
        if makebarplot:
            barplot = bargraph(epitopeproportion, nonepitopeproportion,labs, protein[6], yaxis)

        makeplot = False
        if makeplot:
            makeboxplot = boxplot(masterdata,labels, protein[3])


        appendproteindata = True
        if appendproteindata:
            masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope = appenddata2(epitopeproportion,nonepitopeproportion,masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope)
    labels = ['M1', 'NP']
    outfile = '%s/plots/dnds/usertreecombined_FEL_proportiondndsabove1.pdf' % os.getcwd()
   # print masterproteinproportion_epitope
   # print masterproteinproportion_nonepitope
    M1_epitope = []
    M1_nonepitope = []
    NP_epitope = []
    NP_nonepitope = []
    M1_epitope.append(masterproteinproportion_humanepitope[0])
    M1_epitope.append(masterproteinproportion_swineepitope[0])
    M1_nonepitope.append(masterproteinproportion_humannonepitope[0])
    M1_nonepitope.append(masterproteinproportion_swinenonepitope[0])

    #NP_nonepitope.append(masterproteinproportion_swinenonepitope[0])
    NP_epitope.append(masterproteinproportion_humanepitope[1])
    NP_epitope.append(masterproteinproportion_swineepitope[1])
    NP_nonepitope.append(masterproteinproportion_humannonepitope[1])
    NP_nonepitope.append(masterproteinproportion_swinenonepitope[1])
    print 'plot'

    print masterproteinproportion_humanepitope
    print masterproteinproportion_swineepitope
    print masterproteinproportion_humannonepitope
    print M1_epitope 
    print M1_nonepitope 
    print NP_epitope 
    print NP_nonepitope 
    print 'plot'
    scatter = scatterplot(M1_epitope,NP_epitope,M1_nonepitope,NP_nonepitope,outfile)


   # plotallproteins = True
    #if plotallproteins:
    #print 'plot'
        #plotdndsgreaterthan1 = bargraph4(masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope,labels, outfile, yaxis)
   # scatter = scatterplot(M1_epitope,NP_epitope,M1_nonepitope,NP_nonepitope,outfile)

 
main()