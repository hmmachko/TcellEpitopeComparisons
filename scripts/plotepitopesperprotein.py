'''
This script summarizes the output of epitopefinder, a program we used to map experimentally determined cd8 epitope 
sequences (from Immune Epitope Database) to influenza proteins  NS1, NS2, M1, M2, PB2, PA, NP. 
This script makes two types of graphs,a graph of how many epitopes are in each protein, and a graph for each 
protein that displays the number of epitopes at each amino acid. 

Functions
-------------
*``opencsvfile`` : obtains the number of epitopes per site in a protein from a file
*``countepitopesperprotein`` : obtains the number of epitopes in a protein from a file
*``makedirectory`` : makes a directory if it doesn't exist
*``bargraph`` : bar graph of number of epitpes in each protein
*``listsfromfile`` : 
*``plotbargraph`` :

Input files
-------------
*``combinedepitopesbysite.csv`` : epitopes per site file
*``combinedepitopeslist.csv`` : list of unique epitopes for a protein

Output files
--------------
*``epitopespersite.pdf`` : plot of number of epitopes for each amino acid in a protein
*``numberepitopesperprotein.pdf`` : summary of total epitopes in each protein


'''
import os
import re
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy as np


def countepitopesperprotein(fname,counts):
    '''This function counts the number of unique epitopes in a protein from the output from epitopefinder. 
    *fname* epitopefinder file where each line contains an entry for an epitope
    *counts* list of epitopes for each protein
    
    '''

    fx = open(fname,'r')
    epitopecount = 0
    with fx as f:
        next(f)
        for lines in fx: 
            entry = lines
            epitopecount += 1
    fx.close()
    print epitopecount
    
    counts.append(epitopecount)

def makedirectory(dirname):
    '''This function checks if a directory *dirname* exists, and if it doesn't, it creates it.

    '''
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

def bargraph(data, labs, plotfile, yaxis):
    '''This function creates a bar graph, in this analysis a bar graph of the number of epitopes in a 
    group of proteins.
    *data* a list of the number of epitopes in each protein
    *labs* protein names, used for labeling the bars on the x axis
    *plotfile* name of pdf to save plot
    *yaxis* y axis label
    
    '''
    label_list = []
    for label in labs:
        protein = re.sub('_', ' ',label)
        label_list.append(protein)

    x = xrange(len(data))
    f = pylab.figure()
    bar_width = 0.7
    ax = f.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.bar(x,data,bar_width,align = 'center',color = [95/256.0,158/256.0,209/256.0]) #color = "DodgerBlue"
    ax.set_xticks(x)
    ax.set_xticklabels(label_list,fontsize = 26)
    ax.set_ylabel(yaxis,fontsize = 28)
    matplotlib.pyplot.yticks(fontsize = 24)
    ax.set_xlim([0-bar_width, len(x)])
    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def listsfromfile(fname):
    '''This function returns a list of amino acid sites, a list of number of epitopes per site, and the length of a protein from a file *fname*
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
    return site,eppersites,len(site)

def plotbargraph(xval,yval,protein, plotfile,xaxislim):  
    '''This function generates a bar plot of the number of epitopes at each site along a protein. 
    *xval* list of x-values, in this case the amino acid sites in a protein
    *yval* list of y-values, in this case the number of epitopes at each site
    *plotfile* name of pdf to save plot
    *xaxislim* value for the highest x value, in this case, the last amino acid position
    '''
    print "plotting"
   
    protein = re.sub('_', ' ',protein)
    
    pylab.bar(xval,yval, color='DimGray', linewidth= 0)
    pylab.ylabel('number of epitopes',size = 30)
    pylab.xlabel('site', size = 30)
    matplotlib.pyplot.xlim(1,xaxislim)
    matplotlib.pyplot.yticks(fontsize = 30)
    matplotlib.pyplot.xticks(fontsize = 30)
    matplotlib.pyplot.ylim(0,9.5)
    matplotlib.pyplot.tight_layout()
    pylab.show()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()

def main():
    epitopetypes = ('cd8',)
    proteins = ('NS1','NS2','M2','PB2','PA','M1','NP')

    epsperprotein = []
    for epitope in epitopetypes:
        outputdir = '%s/plots'% (os.getcwd())
        outputsubdir = '%s/plots/%s'%(os.getcwd(), epitope)
        outputsubsubdir = '%s/plots/%s/epitopespersite'%(os.getcwd(), epitope)
        checkdirectory = True
        if checkdirectory:
            directory = makedirectory(outputdir)
            directory = makedirectory(outputsubdir)
            directory = makedirectory(outputsubsubdir)

        for protein in proteins:
            print protein
            infile = ('%s/human/%s/%scombinedepitopeslist.csv' %(os.getcwd(),protein,epitope))
            infile2 = ('%s/human/%s/%scombinedepitopesbysite.csv' %(os.getcwd(),protein,epitope))
            getepitopespersite = True
            if getepitopespersite:
                sitelist,epspersite,lenprotein = listsfromfile(infile2)

            plotepitopespersite = True
            epitopepersiteoutfile ='%s/plots/%s/epitopespersite/%sepitopespersite.pdf' % (os.getcwd(),epitope,protein)
            if plotepitopespersite:
                ploteppersite = plotbargraph(sitelist,epspersite,protein, epitopepersiteoutfile,lenprotein)
            
            countepitopes = True 
            if countepitopes:
                epitopesperprotein = countepitopesperprotein(infile,epsperprotein)
        
        bar = True
        outfile = ('%s/plots/%s/numberepitopesperprotein.pdf' %(os.getcwd(),epitope))
        yax = 'number of unique epitopes'
        if bar:
            graph = bargraph(epsperprotein, proteins, outfile, yax)

main()
