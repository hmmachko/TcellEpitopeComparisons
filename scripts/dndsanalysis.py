'''
This analysis makes plots from 2 types of dN/dS analysis (FUBAR and FEL) for human and swine M1 and NP. 
The first plot is made for both FUBAR and FEL data and is the proportion of sites that have dN/dS > 1 for epitope sites and nonepitope sites. 
The second plot is the cumulative density plot of dN/dS values, and is done for FUBAR analysis.  

Functions
----------
*``readfiletodictionary`` : gets epitopes per site data
*``separateepitopeandnonepitopevalues`` : gets dN/dS values from FUBAR file
*``dataappend`` : formats data for plotting
*``appenddata2`` : formats data for plotting
*``cumulativedensityplot`` : plots cumulative density plot
*``scatterplot`` : makes scatter plot with proportion of sites with dN/dS >1
*``separateepitopeandnonepitopevalues2`` : gets dN/dS values from FEL file
*``makedirectory`` : checks if directory exists and if it doesn't, the function creates it


Input files
-------------
*``protein_dNdStype_report.csv``: dN/dS report. protein is M1 or NP and dNdStype is FUBAR or FEL, located in each host protein folder
*``cd8combinedepitopesbysite.csv`` : epitopes per site file, in host protein folder
*``protein_FUBAR_cumulativedensity_dnds.pdf`` : cumulative density plot, protein is M1 or NP
*``protein_dNdStype_report_fix.csv`` : adjusted format dN/dS report. protein is M1 or NP and dNdStype is FUBAR or FEL, located in each host protein folder
*``proportionhighdnds.pdf`` : plot with proportion of sites that have dN/dS > 1
'''
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

def separateepitopeandnonepitopevalues(filename,epitopedict,outfile,proteinlength,epitopecount,epitopeproportion,nonepitopeproportion,posterior_positivesel_epitope, posterior_positivesel_nonepitope):
    '''This function is for the FUBAR file and first fixes the return in *filename* so that it is recognized by python and creates a new *outfile*. It then makes two lists, a list
    of dN/dS values for epitope sites and a list for nonepitope sites. From this list, this function calculates the proportion of sites (epitope or nonepitope)
    that have dN/dS values greater than 1. Epitope proportion and nonepitope proportions with dN/dS values greater than 1 are appened to a list and returned.
    The lists with epitope and nonepitope dN/dS values are also returned
    *``filename`` : FUBAR dN/dS data file 
    *``epitopedict`` : dictionary with amino acid as key and value as number of epitopes at that site
    *``outfile`` : name of corrected (for carriage return) FUBAR file
    *``proteinlength`` : length of protein
    *``epitopecount`` : number of sites with epitopes
    *``epitopeproportion`` : proportion of sites with epitopes
    *``nonepitopeproportion`` : proportion of sites without epitopes
    '''   
    fx = open(filename, 'r')
    dnds_epitope = []
    dnds_nonepitope = []
    print 'reading file'
    for lines in fx:
        lines = lines.replace('\r','\n')
        #print lines
    f = open(outfile,'w')
    f.write(lines)
    f.close()
    dnds_epitope = []
    dnds_nonepitope = []
    epgreaterthan1 = 0
    nonepgreaterthan1 = 0
    posterior_positiveselection_epitope = []
    posterior_positiveselection_nonepitope = []
    propposterior_epitope = 0
    propposterior_nonepitope = 0
    print 'checking adjusted file'
    f = open(outfile,'r')
    with f as fx:
        next(fx)
        for lines in fx:
            entry = lines.strip().split(',')
            if epitopedict[int(entry[0])] != 0:
                dnds = float(entry[2])/float(entry[1])
                posterior = float(entry[4])
                posterior_positiveselection_epitope.append(float(entry[4]))
                dnds_epitope.append(dnds)
                if dnds >1:
                    epgreaterthan1 +=1
                if posterior >0.9:
                    print entry[0], entry[4]
                    propposterior_epitope +=1
                   
            else:
                dnds = float(entry[2])/float(entry[1])
                dnds_nonepitope.append(dnds)
                posterior = float(entry[4])
                posterior_positiveselection_epitope.append(float(entry[4]))
                if dnds >1:
                    nonepgreaterthan1 +=1
                if posterior >0.9:
                    print entry[0], entry[4]
                    propposterior_nonepitope +=1
                    
    f.close()
    proportiongreaterthan1epitope = float(epgreaterthan1)/epitopecount
    proportiongreaterthan1nonepitope = float(nonepgreaterthan1)/(proteinlength -epitopecount)
    epitopeproportion.append(proportiongreaterthan1epitope)
    nonepitopeproportion.append(proportiongreaterthan1nonepitope)
    proportionposteriorpositiveepitope = float(propposterior_epitope)/epitopecount
    proportionposteriorpositivenonepitope = float(propposterior_nonepitope)/(proteinlength -epitopecount)
    posterior_positivesel_epitope.append(proportionposteriorpositiveepitope) 
    posterior_positivesel_nonepitope.append(proportionposteriorpositivenonepitope)
    print 'epitope posterior proportion: %s , epitope posterior proportion: %s' % (proportionposteriorpositiveepitope,proportionposteriorpositivenonepitope)
    return dnds_epitope,dnds_nonepitope,epitopeproportion, nonepitopeproportion, posterior_positivesel_epitope, posterior_positivesel_nonepitope

def dataappend(masterdata_dnds,human_dnds_epitope,human_dnds_nonepitope,swine_dnds_epitope,swine_dnds_nonepitope):
    ''' This function appends 4 lists (masterdata_dnds,human_dnds_epitope,human_dnds_nonepitope,swine_dnds_epitope,swine_dnds_nonepitope) 
    to 1 master list and returns the master list of lists.'''
    masterdata_dnds.append(human_dnds_epitope)    
    masterdata_dnds.append(human_dnds_nonepitope)
    masterdata_dnds.append(swine_dnds_epitope)
    masterdata_dnds.append(swine_dnds_nonepitope)
    return masterdata_dnds

def appenddata2(dataa,datab,mastera,masterb,masterc,masterd):
    '''This function takes a list and puts the first and second value in each list in a different list. This is done for two lists 
    (*dataa* and *datab*) and the output lists are *mastera*, *masterb*, *masterc*, *masterd*.
    '''
    mastera.append(dataa[0])
    masterb.append(datab[0])
    masterc.append(dataa[1])
    masterd.append(datab[1])
    return mastera, masterb,masterc,masterd

def cumulativedensityplot(masterdata,outfile,xlabels):
    '''This function makes a cumulative density plot from a list of lists *masterdata*, a specified pdf *outfile*, and x labels *xlabels.
    This plot displays a legend based on x labels
    '''
    
    fig, ax = matplotlib.pyplot.subplots()
    index = 0
    for data in masterdata:
        yaxis = []
        numberofdatapoints = len(data)
        sorted_data = np.sort(data)
        question = np.arange(sorted_data.size)
        #print question

        for i,points in enumerate(sorted_data):
            ypoint = (float(i)/numberofdatapoints) + (float(1)/numberofdatapoints)
            yaxis.append(ypoint)
       # print yaxis
       # print len(yaxis)
       # print len(question)
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

def scatterplot(datapointsa,datapointsb,datapointsc,datapointsd,plotfile,yaxlim,yaxtitle,yaxticksend,yaxincrement):
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
    print 'values for plotting'
    print datapointsa,datapointsb,datapointsc,datapointsd
    x_val = [1,2]
    x_val2 = [3,4]
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.borderpad'] =.5 
    matplotlib.pyplot.plot(x_val, a, '-o',label = 'epitope',color = [95/256.0,158/256.0,209/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, b, '-o',color = [95/256.0,158/256.0,209/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val, c, '-o',label = 'nonepitope',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)
    matplotlib.pyplot.plot(x_val2, e, '-o',color = [255/256.0,128/256.0,14/256.0],markersize=15,linewidth=3)

    #print 'values for plotting'
    #print datapointsa,datapointsb,datapointsc,datapointsd
    matplotlib.pyplot.ylabel(yaxtitle,fontsize = 30)
    matplotlib.pyplot.xlabel('M1                     NP',fontsize = 30)
    
    if 'significant'in plotfile:
        print 'posterior'
        matplotlib.pyplot.legend(loc=0, fontsize = 30)
    elif ('FEL' in plotfile) and ('dnds' in plotfile):
        print 'fel and dnds'
        matplotlib.pyplot.legend(loc=4, fontsize = 30)
    else:
        print 'other'
        matplotlib.pyplot.legend(fontsize = 30)
    xlabs = ['human','swine','human','swine','human','swine']
    xtick = [1,2,3,4]
    matplotlib.pyplot.xticks(xtick,xlabs, fontsize = 30)
    matplotlib.pyplot.yticks(fontsize = 30)   
    matplotlib.pyplot.xlim([0.5, 4.5])
    matplotlib.pyplot.ylim(-0.5,yaxticksend)
    matplotlib.pyplot.yticks(np.arange(0,yaxticksend,yaxincrement))
    matplotlib.pyplot.tight_layout()
    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close() 

def separateepitopeandnonepitopevalues2(filename,epitopedict,outfile,proteinlength,epitopecount,epitopeproportion,nonepitopeproportion,epitopeposselproportion,nonepitopeposselproportion):
    '''This function is for the FEL file and first fixes the return in *filename* so that it is recognized by python and creates a new *outfile*. It then makes two lists, a list
    of dN/dS values for epitope sites and a list for nonepitope sites. From this list, this function calculates the proportion of sites (epitope or nonepitope)
    that have dN/dS values greater than 1. Epitope proportion and nonepitope proportions with dN/dS values greater than 1 are appened to a list and returned.
    *``filename`` : FUBAR dN/dS data file 
    *``epitopedict`` : dictionary with amino acid as key and value as number of epitopes at that site
    *``outfile`` : name of corrected (for carriage return) FUBAR file
    *``proteinlength`` : length of protein
    *``epitopecount`` : number of sites with epitopes
    *``epitopeproportion`` : proportion of sites with epitopes
    *``nonepitopeproportion`` : proportion of sites without epitopes
    '''  
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
    proppossel_epitope=0
    proppossel_nonepitope =0
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
                if float(entry[4]) >0:
                    if float(entry[8]) <0.1:
                        print entry[0], entry[8]
                        proppossel_epitope +=1
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
                if float(entry[4]) >0:
                    if float(entry[8]) <0.1:
                        print entry[0], entry[8]
                        proppossel_nonepitope +=1
                    else:
                        continue

    proportiongreaterthan1epitope = float(epgreaterthan1)/epitopecount
    proportiongreaterthan1nonepitope = float(nonepgreaterthan1)/(proteinlength -epitopecount)
    epitopeproportion.append(proportiongreaterthan1epitope)
    nonepitopeproportion.append(proportiongreaterthan1nonepitope)
    #####finish
    proportiongposselepitope = float(proppossel_epitope)/epitopecount
    proportionposselnonepitope = float(proppossel_nonepitope)/(proteinlength -epitopecount)
    epitopeposselproportion.append(proportiongposselepitope)
    nonepitopeposselproportion.append(proportionposselnonepitope)
    totalep = 0
    totalnonep = 0
    for entry in norm_dnds_epitope:
        totalep+=entry
    for entry in norm_dnds_nonepitope:
        totalnonep+= entry
    print 'FEL epitope: %s nonepitope: %s' % (epitopeposselproportion, nonepitopeposselproportion)
    return epitopeproportion, nonepitopeproportion, epitopeposselproportion,nonepitopeposselproportion

def makedirectory(dirname):
    '''This function checks if a directory exists, and if it doesn't, it creates the directory
    *dirname* directory name 
    '''
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
       
def main():
    dnds_analysis = ('FUBAR','FEL')
    seq_format = ('alignedseq','usertree') 
    checkdirectory = True 
    outputsubdir = '%s/plots/dnds'%(os.getcwd())          
    if checkdirectory:
        directory = makedirectory(outputsubdir)
    for analysistype in dnds_analysis:
        for seq in seq_format:


            humanfiles = (
                '%s/dnds/%s/%s/human_M1_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_NP_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                )
            swinefiles = (
                '%s/dnds/%s/%s/swine_M1_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_NP_%s_report.csv' % (os.getcwd(),analysistype,seq,analysistype),
                )
            epitopefiles = (
                '%s/human/M1/cd8combinedepitopesbysite.csv' % os.getcwd(),
                '%s/human/NP/cd8combinedepitopesbysite.csv' % os.getcwd(),
                )

            dndscumulativedensityoutfiles = (
                '%s/plots/dnds/%sM1_FUBAR_cumulativedensity_dnds.pdf' % (os.getcwd(),seq),
                '%s/plots/dnds/%sNP_FUBAR_cumulativedensity_dnds.pdf' % (os.getcwd(),seq),
                )
            adjustedhumanfiles = (    
                '%s/dnds/%s/%s/human_M1_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/human_NP_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                )
            adjustedswinefiles = (      
                '%s/dnds/%s/%s/swine_M1_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                '%s/dnds/%s/%s/swine_NP_%s_report_fix.csv' % (os.getcwd(),analysistype,seq,analysistype),
                )

            proportiondndsgreater1plotoutfiles = (
                '%s/plots/dnds/%s_%s_proportionhighdnds.pdf' % (os.getcwd(),analysistype,seq),

                )
            proportionposteriorsignificantFUBAR = (
                '%s/plots/dnds/%s_%s_proportionposteriorsignificant.pdf' % (os.getcwd(),analysistype,seq),

                )
            proportionsignificantFEL = (
                '%s/plots/dnds/%s_%s_proportionsignificant.pdf' % (os.getcwd(),analysistype,seq),

                )
            masterproteinproportion_humanepitope = []
            masterproteinproportion_humannonepitope = []
            masterproteinproportion_swineepitope = []
            masterproteinproportion_swinenonepitope = []
            M1_masterdata_dnds = []
            NP_masterdata_dnds = []
            M1epitopeproportion = []
            M1nonepitopeproportion = []
            NPepitopeproportion = []
            NPnonepitopeproportion = []
            M1epitopeposteriorproportion = []
            M1nonepitopeposteriorproportion = []
            NPepitopeposteriorproportion = [] 
            NPnonepitopeposteriorproportion = []
            M1posselepitope = [] 
            M1posselnonepitope = []
            NPposselepitope = []
            NPposselnonepitope=[]
            readepitopes = True
            if readepitopes:             
                M1epitopes,M1proteinlength,M1epitopesites = readfiletodictionary(epitopefiles[0], makelist = False)
                NPepitopes,NPproteinlength,NPepitopesites = readfiletodictionary(epitopefiles[1], makelist = False)
            getdndsvalues = True
            if getdndsvalues:
                if analysistype == 'FUBAR':
                    #print humanfiles[0]
                    print adjustedhumanfiles[0]
                    M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion, M1epitopeposteriorproportion, M1nonepitopeposteriorproportion = separateepitopeandnonepitopevalues(humanfiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion,M1epitopeposteriorproportion, M1nonepitopeposteriorproportion)
                    M1_swine_dnds_epitope,M1_swine_dnds_nonepitope,M1epitopeproportion,M1nonepitopeproportion,M1epitopeposteriorproportion, M1nonepitopeposteriorproportion = separateepitopeandnonepitopevalues(swinefiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion,M1epitopeposteriorproportion, M1nonepitopeposteriorproportion)
                    NP_human_dnds_epitope,NP_human_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion,NPepitopeposteriorproportion, NPnonepitopeposteriorproportion = separateepitopeandnonepitopevalues(humanfiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion,NPepitopeposteriorproportion, NPnonepitopeposteriorproportion)
                    NP_swine_dnds_epitope,NP_swine_dnds_nonepitope,NPepitopeproportion,NPnonepitopeproportion,NPepitopeposteriorproportion, NPnonepitopeposteriorproportion = separateepitopeandnonepitopevalues(swinefiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion,NPepitopeposteriorproportion, NPnonepitopeposteriorproportion)
                  
                if analysistype == 'FEL':
                    M1epitopeproportion, M1nonepitopeproportion, M1posselepitope, M1posselnonepitope = separateepitopeandnonepitopevalues2(humanfiles[0],M1epitopes,adjustedhumanfiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion,M1posselepitope, M1posselnonepitope)
                    M1epitopeproportion, M1nonepitopeproportion, M1posselepitope, M1posselnonepitope = separateepitopeandnonepitopevalues2(swinefiles[0],M1epitopes,adjustedswinefiles[0],M1proteinlength,M1epitopesites,M1epitopeproportion,M1nonepitopeproportion,M1posselepitope, M1posselnonepitope)
                    NPepitopeproportion, NPnonepitopeproportion, NPposselepitope, NPposselnonepitope = separateepitopeandnonepitopevalues2(humanfiles[1],NPepitopes,adjustedhumanfiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion, NPposselepitope, NPposselnonepitope)
                    NPepitopeproportion, NPnonepitopeproportion, NPposselepitope, NPposselnonepitope = separateepitopeandnonepitopevalues2(swinefiles[1],NPepitopes,adjustedswinefiles[1],NPproteinlength,NPepitopesites,NPepitopeproportion,NPnonepitopeproportion, NPposselepitope, NPposselnonepitope)
                   
            if analysistype == 'FUBAR':
                appenddata = True
            if appenddata:
                M1_masterdata_dnds = dataappend(M1_masterdata_dnds,M1_human_dnds_epitope,M1_human_dnds_nonepitope,M1_swine_dnds_epitope,M1_swine_dnds_nonepitope)
                NP_masterdata_dnds = dataappend(NP_masterdata_dnds,NP_human_dnds_epitope,NP_human_dnds_nonepitope,NP_swine_dnds_epitope,NP_swine_dnds_nonepitope)
               
            if analysistype =='FUBAR':
                makecumulativedensityplot = True
            xlabels = ['human epitope', 'human nonepitope', 'swine epitope', 'swine nonepitope']
            if makecumulativedensityplot:

                M1densityplot = cumulativedensityplot(M1_masterdata_dnds, dndscumulativedensityoutfiles[0],xlabels)
                NPdensityplot = cumulativedensityplot(NP_masterdata_dnds, dndscumulativedensityoutfiles[1],xlabels)
               
            plotallproteins = True
            if plotallproteins:
                yaxtitle = '% sites with dN/dS > 1'
                yaxlim=[-1,14]
                yaxticksend=13
                yaxincrement=4
                
                scatterplot(M1epitopeproportion,NPepitopeproportion,M1nonepitopeproportion,NPnonepitopeproportion, proportiondndsgreater1plotoutfiles[0],yaxlim,yaxtitle,yaxticksend,yaxincrement)
            if analysistype =='FUBAR':
                plotproportionposterior = True
                yaxlim=[-.5,4]
                yaxticksend=3.5
                yaxincrement=1
                yaxtitle = '% sites under selection'
                if plotproportionposterior:
                    scatterplot(M1epitopeposteriorproportion,NPepitopeposteriorproportion,M1nonepitopeposteriorproportion,NPnonepitopeposteriorproportion, proportionposteriorsignificantFUBAR[0],yaxlim,yaxtitle,yaxticksend,yaxincrement)

            if analysistype =='FEL':
                plotproportionpossel = True
                yaxlim=[-.5,4]
                yaxticksend=3.5
                yaxincrement=1
                yaxtitle = '% sites under selection'
                if plotproportionpossel:


                    scatterplot(M1posselepitope,NPposselepitope,M1posselnonepitope,NPposselnonepitope, proportionsignificantFEL[0],yaxlim,yaxtitle,yaxticksend,yaxincrement)   


main()









