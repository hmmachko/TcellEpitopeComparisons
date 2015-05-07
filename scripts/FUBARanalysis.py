
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

    posterior_positiveselection_epitope = [] 
    posterior_positiveselection_nonepitope = []
    dnds_epitope = []
    dnds_nonepitope = []
    epgreaterthan1 = 0
    nonepgreaterthan1 = 0

    print 'checking adjusted file'
    f = open(outfile,'r')
    with f as fx:
        next(fx)
        for lines in fx:
            entry = lines.strip().split(',')
            #print entry
            #print entry[3]
            if epitopedict[int(entry[0])] != 0:
                posterior_positiveselection_epitope.append(float(entry[4]))
                dnds = float(entry[2])/float(entry[1])
                #print dnds
                dnds_epitope.append(dnds)
                if dnds >1:
                    epgreaterthan1 +=1
            else:
                posterior_positiveselection_nonepitope.append(float(entry[4]))
                dnds = float(entry[2])/float(entry[1])
                #print dnds
                dnds_nonepitope.append(dnds)
                if dnds >1:
                    nonepgreaterthan1 +=1
    f.close()

    proportiongreaterthan1epitope = float(epgreaterthan1)/epitopecount
    proportiongreaterthan1nonepitope = float(nonepgreaterthan1)/(proteinlength -epitopecount)
    epitopeproportion.append(proportiongreaterthan1epitope)
    nonepitopeproportion.append(proportiongreaterthan1nonepitope)

    return posterior_positiveselection_epitope,posterior_positiveselection_nonepitope, dnds_epitope,dnds_nonepitope,epitopeproportion, nonepitopeproportion

def makeviolinplot(datalists,plotfile,xlabels,ylab):
    '''This function constructs a violin plot using matplotlib

    *datalists* lists of datapoints
    *plotfile* outputfile of violinplot in pdf format
    *protein* xlabel
    *masterdeltaf* list of points, 1 placed on each violin plot 

    '''
    pos =  [1,2,4,5]
    fig, ax = matplotlib.pyplot.subplots()
    plot = matplotlib.pyplot.violinplot(datalists, pos, points=60, widths=0.4, showmeans=True,
                      showextrema=False, showmedians=False, bw_method="scott")
    #plot = matplotlib.pyplot.scatter(pos,masterdeltaf,s=60, color='Red')
    ax.set_xticks(pos)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel(ylab)
    matplotlib.pyplot.ylim([-0.1, 2])
    #ax.gca.set_ylabel(r'$\boldsymbol{\delta}average epitope change per mutation$')

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def dataappend(masterdata_dnds,human_dnds_epitope,human_dnds_nonepitope,swine_dnds_epitope,swine_dnds_nonepitope):
    masterdata_dnds.append(human_dnds_epitope)    
    masterdata_dnds.append(human_dnds_nonepitope)
    masterdata_dnds.append(swine_dnds_epitope)
    masterdata_dnds.append(swine_dnds_nonepitope)



    return masterdata_dnds

def appenddata2(dataa,datab,mastera,masterb,masterc,masterd):
   

    mastera.append(dataa[0])
    masterb.append(datab[0])
    masterc.append(dataa[1])
    masterd.append(datab[1])
    

    return mastera, masterb,masterc,masterd

def bargraph(data, data2,labs, plotfile, yaxis):
   
    x = xrange(len(data))

    n_groups = len(data)
    index = np.arange(n_groups)
    bar_width = 0.35

    fig, ax = matplotlib.pyplot.subplots()
    aplot = matplotlib.pyplot.bar(index,data, bar_width, label = 'epitope')#color = "DodgerBlue"
    bplot = matplotlib.pyplot.bar(index+ bar_width,data2, bar_width, color = 'Red', label = 'nonepitope')#color = "FireBrick"
    #matplotlib.pyplot.xlabel('protein')
    matplotlib.pyplot.ylabel(yaxis,fontsize = 15)
    #matplotlib.pyplot.title(epitopetype)
    matplotlib.pyplot.xticks(index + bar_width, labs,fontsize = 15)
    matplotlib.pyplot.yticks(fontsize = 11)
    matplotlib.pyplot.legend(fontsize = 13)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

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
    matplotlib.pyplot.legend(fontsize = 13)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def cumulativedensityplot(masterdata,outfile,xlabels):
    #dataarray = np.array(data)
    fig, ax = matplotlib.pyplot.subplots()
    index = 0
    for data in masterdata:
        yaxis = []
        numberofdatapoints = len(data)
        sorted_data = np.sort(data)
        question = np.arange(sorted_data.size)
        print question

        for i,points in enumerate(sorted_data):
            ypoint = (float(i)/numberofdatapoints) + (float(1)/numberofdatapoints)
            yaxis.append(ypoint)
        print yaxis
        print len(yaxis)
        print len(question)
        #print sorted_data
        plota = matplotlib.pyplot.step(sorted_data, yaxis,label = xlabels[index])
        index +=1
        #np.arange(sorted_data.size)
    matplotlib.pyplot.ylabel('cumulative percent',fontsize = 15)
    matplotlib.pyplot.yticks(fontsize = 11)
    matplotlib.pyplot.xticks(fontsize = 11)
    matplotlib.pyplot.xlabel('dnds',fontsize = 15)

    matplotlib.pyplot.legend(fontsize = 13)

    pylab.show()
    pylab.savefig(outfile)
    pylab.clf()
    pylab.close()

def scatterplot(datapointsa,datapointsb,datapointsc,datapointsd,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk


    #x_val = [1,2,4,5]
    x_val = [1,2]
    x_val2 = [3,4]
    #fig = matplotlib.pyplot.figure()
   # ax1 = fig.add_subplot(111)
    matplotlib.pyplot.plot(x_val, datapointsa,'-o',ms = 10,color = [95/256.0,158/256.0,209/256.0], label = 'epitope')
    matplotlib.pyplot.plot(x_val2, datapointsb,'-o',ms = 10,color = [95/256.0,158/256.0,209/256.0])
    matplotlib.pyplot.plot(x_val, datapointsc,'-o',ms = 10,color =[255/256.0,128/256.0,14/256.0],label = 'nonepitope')
    matplotlib.pyplot.plot(x_val2, datapointsd,'-o',ms = 10, color = [255/256.0,128/256.0,14/256.0])


    #matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'epitope',s = 70)
   # matplotlib.pyplot.scatter(x_val2, datapointsb,color = [95/256.0,158/256.0,209/256.0],s = 70)
   # matplotlib.pyplot.scatter(x_val, datapointsc,color = [255/256.0,128/256.0,14/256.0],label = 'nonepitope',s = 70,marker = '*')
   # matplotlib.pyplot.scatter(x_val2, datapointsd,color = [255/256.0,128/256.0,14/256.0],s = 70,marker = '*')
    matplotlib.pyplot.ylabel('fraction of sites with dnds > 1',fontsize = 17)
    matplotlib.pyplot.xlabel('M1                                     NP',fontsize = 17)
    matplotlib.pyplot.legend(loc=1, fontsize = 15)
    xlabs = ['human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 17)
    matplotlib.pyplot.yticks(fontsize = 13)
    matplotlib.pyplot.ylim([-0.01, 0.07])
    matplotlib.pyplot.xlim([0.5, 4.5])

    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close() 

def main():

    humanFUBARfiles = (
        #'%s/dnds/FUBAR/human_PA_FUBAR_report.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/human_M1_FUBAR_report.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/human_NP_FUBAR_report.csv' % os.getcwd(),
        )
    swineFUBARfiles = (
       # '%s/dnds/FUBAR/swine_PA_FUBAR_report.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/swine_M1_FUBAR_report.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/swine_NP_FUBAR_report.csv' % os.getcwd(),
        )
    epitopefiles = (
      #  '%s/human/PA/cd8combinedepitopesbysite.csv' % os.getcwd(),
        '%s/human/M1/cd8combinedepitopesbysite.csv' % os.getcwd(),
        '%s/human/NP/cd8combinedepitopesbysite.csv' % os.getcwd(),
        )
    posteriorviolinplotoutfiles = (
      #  '%s/plots/dnds/PA_FUBAR_posterior_positiveselection.pdf' % os.getcwd(),
        '%s/plots/dnds/NP_FUBAR_posterior_positiveselection.pdf' % os.getcwd(),
        '%s/plots/dnds/NP_FUBAR_dndsviolin.pdf' % os.getcwd(),
        )
    dndscumulativedensityoutfiles = (
      #  '%s/plots/dnds/PA_FUBAR_cumulativedensity_dnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertreeM1_FUBAR_cumulativedensity_dnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertreeNP_FUBAR_cumulativedensity_dnds.pdf' % os.getcwd(),
        )
    adjustedhumanFUBARfiles = (
      #  '%s/dnds/FUBAR/human_PA_FUBAR_report_fix.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/human_M1_FUBAR_report_fix.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/human_NP_FUBAR_report_fix.csv' % os.getcwd(),
        )
    adjustedswineFUBARfiles = (
       # '%s/dnds/FUBAR/swine_PA_FUBAR_report_fix.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/swine_M1_FUBAR_report_fix.csv' % os.getcwd(),
        '%s/dnds/FUBAR/usertree/swine_NP_FUBAR_report_fix.csv' % os.getcwd(),
        )

    proportiondndsgreater1plotoutfiles = (
       # '%s/plots/dnds/PA_FUBAR_proportionhighdnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertreeM1_FUBAR_proportionhighdnds.pdf' % os.getcwd(),
        '%s/plots/dnds/usertreeNP_FUBAR_proportionhighdnds.pdf' % os.getcwd(),
        )

    masterproteinproportion_humanepitope = []
    masterproteinproportion_humannonepitope = []
    masterproteinproportion_swineepitope = []
    masterproteinproportion_swinenonepitope = []

    #analysis = list(zip(humanFUBARfiles,swineFUBARfiles,epitopefiles,posteriorviolinplotoutfiles,dndscumulativedensityoutfiles,adjustedhumanFUBARfiles,adjustedswineFUBARfiles,proportiondndsgreater1plotoutfiles))
    #for protein in analysis:
    M1_masterdata_dnds = []
    NP_masterdata_dnds = []
    masterdata_posterior = []
    M1epitopeproportion = []
    M1nonepitopeproportion = []
    NPepitopeproportion = []
    NPnonepitopeproportion = []
    readepitopes = True
    if readepitopes:
        M1epitopes,M1proteinlength,M1epitopesites = readfiletodictionary(epitopefiles[0], makelist = False)
        NPepitopes,NPproteinlength,NPepitopesites = readfiletodictionary(epitopefiles[1], makelist = False)


    getdndsvalues = True
    if getdndsvalues:
        M1_human_posterior_epitope,M1_human_posterior_nonepitope,M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion = separateepitopeandnonepitopevalues(humanFUBARfiles[0],M1epitopes,adjustedhumanFUBARfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
        M1_swine_posterior_epitope,M1_swine_posterior_nonepitope,M1_swine_dnds_epitope,M1_swine_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion = separateepitopeandnonepitopevalues(swineFUBARfiles[0],M1epitopes,adjustedhumanFUBARfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
        NP_human_posterior_epitope,NP_human_posterior_nonepitope,NP_human_dnds_epitope,NP_human_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion = separateepitopeandnonepitopevalues(humanFUBARfiles[1],NPepitopes,adjustedhumanFUBARfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)
        NP_swine_posterior_epitope,NP_swine_posterior_nonepitope,NP_swine_dnds_epitope,NP_swine_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion = separateepitopeandnonepitopevalues(swineFUBARfiles[1],NPepitopes,adjustedhumanFUBARfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)

    appenddata = True
    if appenddata:
        M1_masterdata_dnds = dataappend(M1_masterdata_dnds,M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1_swine_dnds_epitope,M1_swine_dnds_nonepitope)
        NP_masterdata_dnds = dataappend(NP_masterdata_dnds,NP_human_dnds_epitope,NP_human_dnds_nonepitope,NP_swine_dnds_epitope,NP_swine_dnds_nonepitope)


    violinplot_posterior = True
    ylab = 'posterior probability for positive selection'
    xlabels = ['human epitope', 'human nonepitope', 'swine epitope', 'swine nonepitope']
    if violinplot_posterior:
        masterdata_posterior.append(NP_human_posterior_epitope)
        masterdata_posterior.append(NP_human_posterior_nonepitope)
        masterdata_posterior.append(NP_swine_posterior_epitope)
        masterdata_posterior.append(NP_swine_posterior_nonepitope)

        vplot = makeviolinplot(masterdata_posterior, posteriorviolinplotoutfiles[0],xlabels,ylab)

    violinplot_dnds = True
    ylab = 'dnds'
    if violinplot_dnds:
        #master_dnds = []
        #masterdata_dnds.append 
        vplot = makeviolinplot(NP_masterdata_dnds, posteriorviolinplotoutfiles[1],xlabels,ylab)

    makebarplot = False
    labs = ['human', 'swine']
    yaxis = ('fraction of sites with dnds > 1')
    if makebarplot:
        barplot = bargraph(epitopeproportion, nonepitopeproportion,labs, protein[7], yaxis)

    makecumulativedensityplot = True
    if makecumulativedensityplot:
        M1densityplot = cumulativedensityplot(M1_masterdata_dnds, dndscumulativedensityoutfiles[0],xlabels)
        NPdensityplot = cumulativedensityplot(NP_masterdata_dnds, dndscumulativedensityoutfiles[1],xlabels)
    


    appendproteindata = False
    if appendproteindata:
        masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope = appenddata2(epitopeproportion,nonepitopeproportion,masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope)
    labels = ['M1', 'NP']
    outfile = '%s/plots/dnds/usertreecombined_FUBAR_proportiondndsabove1.pdf' % os.getcwd()
   # print masterproteinproportion_epitope
   # print masterproteinproportion_nonepitope

    plotallproteins = True
    if plotallproteins:
        #plotdndsgreaterthan1 = bargraph4(masterproteinproportion_humanepitope,masterproteinproportion_humannonepitope,masterproteinproportion_swineepitope,masterproteinproportion_swinenonepitope,labels, outfile, yaxis)
        scatterplot(M1epitopeproportion,NPepitopeproportion,M1nonepitopeproportion,NPnonepitopeproportion,outfile)


main()









