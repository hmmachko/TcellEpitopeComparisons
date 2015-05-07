import os
import re
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np

def opencsvfile(fname):
 
    proteins = []
    epitopenumber = []
    fx = open(fname,'r')
    with fx as f:
        #next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            epitopenumber.append(int(entry[1]))
            proteins.append(entry[0])
    return proteins, epitopenumber

def countepitopesperprotein(fname,counts):
    fx = open(fname,'r')
    epitopecount = 0
    with fx as f:
        next(f)
        for lines in fx: 
            entry = lines
            epitopecount += 1
    fx.close()
    counts.append(epitopecount)

def makedirectory(dirname):
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

def bargraph(data, labs, plotfile, yaxis,epitopetype):
    x = xrange(len(data))
    f = pylab.figure()
    bar_width = 0.7
    ax = f.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.bar(x,data,bar_width,align = 'center',color = "DodgerBlue")
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylabel(yaxis)
    ax.set_xlabel('protein')
    ax.set_title(epitopetype)
    ax.set_xlim([0-bar_width, len(x)])

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def listsfromfile(fname):
    '''This function returns a list of aminoacid sites and number of epitopes per site from a file *fname*
    *fname* file that contains a titleline as the first line and subsequent lines with sites as first entry and epitopes
    per site as second entry
    ''' 
    site = []
    eppersites = []
    fx = open(fname,'r')
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            site.append(int(entry[0]))
            eppersites.append(float(entry[1]))
    return site,eppersites

def plotbargraph(xval,yval,protein, plotfile):  ####fix
    print "plotting"
   
    protein = re.sub('_', ' ',protein)
    
    pylab.bar(xval,yval, color='DimGray', linewidth= 0)
    pylab.ylabel('number epitopes')
    pylab.xlabel('site')
    pylab.title(protein)

    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def main():
    epitopetypes = ('cd8',)
    proteins = ('NA','NS1','NS2','HA','M2','PB2','PA','M1','PB1','NP')
    epitopesperprotein = []
    for epitope in epitopetypes:
        for protein in proteins:
            infile = ('%s/human/%s/%scombinedepitopeslist.csv' %(os.getcwd(),protein,epitope))
            countepitopes = True
            if countepitopes:
                epitopesperprotein = countepitopesperprotein(infile,epitopesperprotein)

            infile2 = ('%s/human/%s/%scombinedepitopesbysite.csv' %(os.getcwd(),protein,epitope))
            getepitopespersite = True
            if getepitopespersite:
                sitelist,epspersite = listsfromfile(infile2)

            plotepitopespersite = True
            epitopepersiteoutfile ='%s/plots/epitopespersite/%s%sepitopespersite.pdf' % (os.getcwd(),protein,epitope))
            if plotepitopespersite:
                ploteppersite = plotbargraph(sitelist,epspersite,protein, epitopepersiteoutfile)

        checkdirectory = True
        dirname = ('%s/plots/%s'%(os.getcwd(),epitope))
        if checkdirectory:
            directory = makedirectory()
        bar = True
        outfile = ('%s/plots/%s/numberepitopesperprotein.pdf' %(os.getcwd(),epitope))
        yax = 'number of unique epitopes'
        if bar:
            graph = bargraph(epitopesperprotein, proteins, outfile, yax,epitope)

main()
