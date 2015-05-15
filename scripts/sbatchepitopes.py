'''This script submits to sbatch the program epitope finder to align CD8 T-cell epitopes to influenza proteins.
This analysis is done for NS1, NS2, PB1,  M1, M2, PA, and NP. Both human and swine influenza sequences are included.

Functions
-----------
*``makeinfile`` : creates epitope finder infile
*``CreateSbatchFile`` : creates sbatch file

Input files
-------------
*``prot_aligned.fasta`` : human and swine influenza sequences used to align epitopes, a file exists for each protein
*``IEDB_Influenza_Tcell_compact_2014-10-09.csv``: influenza T-cell epitopes frome Immune Epitope Database
*``supertype_classification.txt`` : MHC classification text

Output files
--------------
*``epitopefinderinfile.txt`` : input file with required parameters to run epitope finder
*``runepitopefinder.sbatch`` : sbatch file to run epitope finder
*``combinedepitopeslist.csv`` : file containing epitope entries from ``IEDB_Influenza_Tcell_compact_2014-10-09.csv`` that align to a protein
*``combinedepitopesbysite.csv`` : file containing the number of unique epitopes that map to each amino acid site of a protein. 

Written by Heather Machkovech Feb 2015

'''

import os
import sys
import listFASTA

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def makeinfile(folder, mhcdictionary,infilename,epitopetype):
    '''This function makes the input file to run epitope finder. 
    *``folder`` : folder to write output of epitope finder
    *``mhcdictionary`` : dictionary that contains input parameters for epitope finder
    *``infilename`` : name of input file for epitope finder
    *``epitopetype`` : epitope type
    '''
    fx = open(infilename, 'w')

    fx.write('# input file for epitopefinder_getepitopes.py\n')
    fx.write('iedbfile %s/inputfiles/IEDB_Influenza_Tcell_compact_2014-10-09.csv\n' % os.getcwd())
    fx.write('supertypesfile %s/inputfiles/supertype_classification.txt\n' % os.getcwd())
    fx.write('mhcclass %s\n' % (mhcdictionary['%s'% epitopetype]['mhcclass'])) 
    fx.write('epitopelength %s %s\n'% (mhcdictionary['%s'% epitopetype]['length'][0],mhcdictionary['%s'% epitopetype]['length'][1]))
    fx.write('targetprotsfile %scombined_prot_aligned.fasta\n' % folder)
    fx.write('maxmismatches 1\n')
    fx.write('musclepath /usr/bin/\n')
    fx.write('purgeredundant MHCgroup\n')
    fx.write('purgeredundantoverlap %s\n' % (mhcdictionary['%s'% epitopetype]['purgeredundantoverlap']))
    fx.write('epitopeslistfile %s%scombinedepitopeslist.csv\n' % (folder, epitopetype))
    fx.write('epitopesbysitefile %s%scombinedepitopesbysite.csv\n' % (folder,epitopetype))
    fx.close()

def CreateSbatchFile(fname,time,command):
    '''This function creates an sbatch file 
    *``fname`` : sbatch.run filename
    *``time`` : run time for sbatch in format 000 (hours)
    *``command`` : command for sbatch in format: 'program infile'
    '''
    f = open(fname,'w')
    f.write('#!/bin/sh\n')
    f.write('#SBATCH\n')
    f.write('#PBS -l walltime=%s:00:00\n' % time)
    f.write('%s'%command)
    f.close()
    

def main():

    mhcdict = AutoVivification()
    mhcdict['cd8']['mhcclass'] = 'I'
    mhcdict['cd8']['length'] = [8, 12]
    mhcdict['cd8']['purgeredundantoverlap'] = 8
    mhcdict['cd4']['mhcclass'] = 'II'
    mhcdict['cd4']['length'] = [13, 24]
    mhcdict['cd4']['purgeredundantoverlap'] = 13

    epitopetypes = ('cd8',)

    humanfolders = (
        '%s/human/NS1/prot_aligned.fasta' % os.getcwd(),
        '%s/human/NS2/prot_aligned.fasta' % os.getcwd(),
        '%s/human/M1/prot_aligned.fasta' % os.getcwd(),
        '%s/human/M2/prot_aligned.fasta' % os.getcwd(),
        '%s/human/NP/prot_aligned.fasta' % os.getcwd(),       
        '%s/human/PA/prot_aligned.fasta' % os.getcwd(),       
        '%s/human/PB2/prot_aligned.fasta' % os.getcwd()
        )
    swinefolders = (
        '%s/swine/NS1/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/NS2/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/M1/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/M2/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/NP/prot_aligned.fasta' % os.getcwd(),     
        '%s/swine/PA/prot_aligned.fasta' % os.getcwd(),   
        '%s/swine/PB2/prot_aligned.fasta' % os.getcwd()
        )

    human_swine_seq = list(zip(humanfolders,swinefolders))
    for epitopetype in epitopetypes:
        for pair in human_swine_seq:
            base = pair[0][:-18]
            joinedseq = '%scombined_prot_aligned.fasta' % base
            concatseqs = True
            if concatseqs:
                concatseq = os.system('cat %s %s > %s' % (pair[0], pair[1],joinedseq))


            home = os.path.expanduser("~")
            epitopefile = '%s%sepitopefinderinfile.txt' %(base,epitopetype)
            program = 'epitopefinder_getepitopes.py'
            command = program + ' ' + epitopefile 
            time = 024
            sbatchfname = '%srunepitopefinder.sbatch' % base
            
            make_epinfile = True
            if make_epinfile:
                inputfile = makeinfile(base, mhcdict,epitopefile,epitopetype)

            make_sbatch = True
            if make_sbatch:
                sbatchfile = CreateSbatchFile(sbatchfname,time,command)
    
            run_sbatch = True
            if run_sbatch:
                os.system('sbatch %s' %(sbatchfname)) 

main()