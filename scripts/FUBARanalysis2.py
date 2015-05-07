
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
    matplotlib.pyplot.ylabel('cumulative percent',fontsize = 20)
    matplotlib.pyplot.yticks(fontsize = 18)
    matplotlib.pyplot.xticks(fontsize = 18)
    matplotlib.pyplot.xlabel('dnds',fontsize = 20)

    matplotlib.pyplot.legend(fontsize = 20)

    pylab.show()
    pylab.savefig(outfile)
    pylab.clf()
    pylab.close()

def scatterplot(datapointsa,datapointsb,datapointse,datapointsc,datapointsd,datapointsf,plotfile):
    #fhumantree,fhumantrunk,fswinetree,fswinetrunk
    x=0
    xo=3
    y=0
    yo=4
    x,xo,y,yo = matplotlib.pyplot.axis()
    plot_margin=.25
    a = []
    b = []
    c= []
    e = []
    for d in datapointsa:
        d = 100*d
        a.append(d)
    for d in datapointsb:
        d = 100*d
        b.append(d)
    for d in datapointsc:
        d = 100*d
        c.append(d)
    for d in datapointsd:
        d = 100*d
        e.append(d)

    

    #x_val = [1,2,4,5]
    x_val = [1,2]
    x_val2 = [3,4]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 

    #x_val3 = [5,6]
    matplotlib.pyplot.plot(x_val, a, '-o',label = 'epitope',color = [95/256.0,158/256.0,209/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, b, '-o',color = [95/256.0,158/256.0,209/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val, c, '-o',label = 'nonepitope',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, e, '-o',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)
    #matplotlib.pyplot.plot(x_val3, datapointse, '-o',color = [95/256.0,158/256.0,209/256.0])
    #matplotlib.pyplot.plot(x_val3, datapointsf, '-o',color = [255/256.0,128/256.0,14/256.0])
    print 'values for plotting'
    #matplotlib.pyplot.axis((0 - plot_margin,
         # 3 + plot_margin,
         # 0 - plot_margin,
         # 4 + plot_margin))
    print datapointsa,datapointsb,datapointse,datapointsc,datapointsd,datapointsf
    #matplotlib.pyplot.scatter(x_val, datapointsa,color = [95/256.0,158/256.0,209/256.0], label = 'tree')
    #matplotlib.pyplot.scatter(x_val, datapointsb,color = [255/256.0,128/256.0,14/256.0], label = 'trunk')
    matplotlib.pyplot.ylabel('% sites with dN/dS > 1',fontsize = 30)
    matplotlib.pyplot.xlabel('M1                     NP',fontsize = 30)
    if 'FEL' in plotfile:

        matplotlib.pyplot.legend(loc=4, fontsize = 30)
    else:
        matplotlib.pyplot.legend(fontsize = 30)
    xlabs = ['human','swine','human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 30)
    
    
    matplotlib.pyplot.xlim([0.5, 4.5])
    matplotlib.pyplot.ylim(-1,14)
    matplotlib.pyplot.yticks(np.arange(0, 13, 4))
    matplotlib.pyplot.tight_layout()

    #matplotlib.pyplot.xticks(index + 2*bar_width, xlabs, fontsize = 15)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close() 

def separateepitopeandnonepitopevalues2(filename,epitopedict,outfile,proteinlength,epitopecount,epitopeproportion,nonepitopeproportion):
    
    fx = open(filename, 'r')
    dnds_epitope = []
    dnds_nonepitope = []
    norm_dnds_epitope = []
    norm_dnds_nonepitope = []
    print 'reading file'
    for lines in fx:
        lines = lines.replace('\r','\n')

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

            else:
                if entry[3] != 'Undefined':

                    if entry[3] == 'Infinite':
                        nonepgreaterthan1 +=1
                    elif float(entry[3]) > 1:
                        nonepgreaterthan1 +=1
                    else:
                        continue
    print 'partway there'
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

    return epitopeproportion, nonepitopeproportion

def main():
    dnds_analysis = ('FUBAR','FEL')
    seq_format = ('alignedseq','usertree')  #'usertree'
    for analysistype in dnds_analysis:
        for seq in seq_format:


            humanfiles = (
                '%s/dnds/%s/%s/human_M1_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_NP_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_HA_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype)
                )
            swinefiles = (
                '%s/dnds/%s/%s/swine_M1_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_NP_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_HA_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype)
                )
            epitopefiles = (
                '%s/human/M1/cd8combinedepitopesbysite.csv' % os.getcwd(),
                '%s/human/NP/cd8combinedepitopesbysite.csv' % os.getcwd(),
                '%s/human/HA_H3/antibodyepitopesbysite.csv' % os.getcwd()
                )

            dndscumulativedensityoutfiles = (
                '%s/plots/dnds/%sM1_FUBAR_cumulativedensity_dnds.pdf' % (os.getcwd(),seq),
                '%s/plots/dnds/%sNP_FUBAR_cumulativedensity_dnds.pdf' % (os.getcwd(),seq),
                '%s/plots/dnds/%sHA_FUBAR_cumulativedensity_dnds.pdf' % (os.getcwd(),seq)
                )
            adjustedhumanfiles = (    
                '%s/dnds/%s/%s/human_M1_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_NP_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_NA_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype)
                )
            adjustedswinefiles = (      
                '%s/dnds/%s/%s/swine_M1_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_NP_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_HA_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype)
                )

            proportiondndsgreater1plotoutfiles = (
                '%s/plots/dnds/%s_%s_proportionhighdnds.pdf' % (os.getcwd(),analysistype,seq),
               # '%s/plots/dnds/%sNP_%s_proportionhighdnds.pdf' % (os.getcwd(),analysistype,seq),
               # '%s/plots/dnds/%sH3_%s_proportionhighdnds.pdf' % (os.getcwd(),analysistype,seq)

                )

            masterproteinproportion_humanepitope = []
            masterproteinproportion_humannonepitope = []
            masterproteinproportion_swineepitope = []
            masterproteinproportion_swinenonepitope = []

    #analysis = list(zip(humanFUBARfiles,swineFUBARfiles,epitopefiles,posteriorviolinplotoutfiles,dndscumulativedensityoutfiles,adjustedhumanFUBARfiles,adjustedswineFUBARfiles,proportiondndsgreater1plotoutfiles))
    #for protein in analysis:
            M1_masterdata_dnds = []
            NP_masterdata_dnds = []
            H3_masterdata_dnds = []
        #masterdata_posterior = []
            M1epitopeproportion = []
            M1nonepitopeproportion = []
            H3epitopeproportion = []
            H3nonepitopeproportion = []
            NPepitopeproportion = []
            NPnonepitopeproportion = []
            readepitopes = True
            if readepitopes:             
                M1epitopes,M1proteinlength,M1epitopesites = readfiletodictionary(epitopefiles[0], makelist = False)
                NPepitopes,NPproteinlength,NPepitopesites = readfiletodictionary(epitopefiles[1], makelist = False)
                H3epitopes,H3proteinlength,H3epitopesites = readfiletodictionary(epitopefiles[2], makelist = False)
            getdndsvalues = True
            if getdndsvalues:
                if analysistype == 'FUBAR':
                    print humanfiles[0]
                    print adjustedhumanfiles[0]
                    M1_human_posterior_epitope,M1_human_posterior_nonepitope,M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion = separateepitopeandnonepitopevalues(humanfiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
                    M1_swine_posterior_epitope,M1_swine_posterior_nonepitope,M1_swine_dnds_epitope,M1_swine_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion = separateepitopeandnonepitopevalues(swinefiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
                    NP_human_posterior_epitope,NP_human_posterior_nonepitope,NP_human_dnds_epitope,NP_human_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion = separateepitopeandnonepitopevalues(humanfiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)
                    NP_swine_posterior_epitope,NP_swine_posterior_nonepitope,NP_swine_dnds_epitope,NP_swine_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion = separateepitopeandnonepitopevalues(swinefiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)
                    H3_human_posterior_epitope,H3_human_posterior_nonepitope,H3_human_dnds_epitope,H3_human_dnds_nonepitope,H3epitopeproportion,H3nonepitopeproportion = separateepitopeandnonepitopevalues(humanfiles[2],H3epitopes,adjustedhumanfiles[2],H3proteinlength,H3epitopesites,H3epitopeproportion,H3nonepitopeproportion)
                    H3_swine_posterior_epitope,H3_swine_posterior_nonepitope,H3_swine_dnds_epitope,H3_swine_dnds_nonepitope,H3epitopeproportion,H3nonepitopeproportion = separateepitopeandnonepitopevalues(swinefiles[2],H3epitopes,adjustedhumanfiles[2],H3proteinlength,H3epitopesites,H3epitopeproportion,H3nonepitopeproportion)
                     
                if analysistype == 'FEL':
                    M1epitopeproportion, M1nonepitopeproportion = separateepitopeandnonepitopevalues2(humanfiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
                    M1epitopeproportion, M1nonepitopeproportion = separateepitopeandnonepitopevalues2(swinefiles[0],M1epitopes,adjustedswinefiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion)
                    NPepitopeproportion, NPnonepitopeproportion = separateepitopeandnonepitopevalues2(humanfiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)
                    NPepitopeproportion, NPnonepitopeproportion = separateepitopeandnonepitopevalues2(swinefiles[1],NPepitopes,adjustedswinefiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion)
                    H3epitopeproportion, H3nonepitopeproportion = separateepitopeandnonepitopevalues2(humanfiles[2],H3epitopes,adjustedhumanfiles[2],H3proteinlength,H3epitopesites,H3epitopeproportion,H3nonepitopeproportion)
                    H3epitopeproportion, H3nonepitopeproportion = separateepitopeandnonepitopevalues2(swinefiles[2],H3epitopes,adjustedswinefiles[2],H3proteinlength,H3epitopesites,H3epitopeproportion,H3nonepitopeproportion)
            if analysistype == 'FUBAR':
                appenddata = True
            if appenddata:
                M1_masterdata_dnds = dataappend(M1_masterdata_dnds,M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1_swine_dnds_epitope,M1_swine_dnds_nonepitope)
                NP_masterdata_dnds = dataappend(NP_masterdata_dnds,NP_human_dnds_epitope,NP_human_dnds_nonepitope,NP_swine_dnds_epitope,NP_swine_dnds_nonepitope)
                H3_masterdata_dnds = dataappend(H3_masterdata_dnds,H3_human_dnds_epitope,H3_human_dnds_nonepitope,H3_swine_dnds_epitope,H3_swine_dnds_nonepitope)

            if analysistype =='FUBAR':
                makecumulativedensityplot = True
            xlabels = ['human epitope', 'human nonepitope', 'swine epitope', 'swine nonepitope']
            if makecumulativedensityplot:

                M1densityplot = cumulativedensityplot(M1_masterdata_dnds, dndscumulativedensityoutfiles[0],xlabels)
                NPdensityplot = cumulativedensityplot(NP_masterdata_dnds, dndscumulativedensityoutfiles[1],xlabels)
                H3densityplot = cumulativedensityplot(H3_masterdata_dnds, dndscumulativedensityoutfiles[2],xlabels)

            plotallproteins = True
            if plotallproteins:
                scatterplot(M1epitopeproportion,NPepitopeproportion,H3epitopeproportion,M1nonepitopeproportion,NPnonepitopeproportion,H3nonepitopeproportion, proportiondndsgreater1plotoutfiles[0])


main()









