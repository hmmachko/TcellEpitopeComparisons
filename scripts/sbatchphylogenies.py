'''This script submits for sbatch the command to build combine human and swine influenza BEAST trees for M1 and NP. 

Functions
-----------
*``CreateSbatchFile`` : creates sbatch file

Input files
-------------
*``beastinfile.xml`` : input file to run BEAST

Output files
--------------
*``run.sbatch`` : sbatch file to run BEAST
*``comb_prot_aligned.trees`` : BEAST output file containing thinned trees 
*``comb_prot_aligned.log` : BEAST log output

Written by Heather Machkovech Feb 2015

'''

import os
import sys

def CreateSbatchFile(fname,time,command):
    '''This function creates an sbatch file with filename *fname*, time for running *time* (in hours) and sbatch command *command*
    '''
    f = open(fname,'w')
    f.write('#!/bin/sh\n')
    f.write('#SBATCH\n')
    f.write('#PBS -l walltime=%s:00:00\n' % time)
    f.write('%s'%command)
    f.close()


def main():

    folders = (
        '%s/human/M1/' % os.getcwd(),
        '%s/human/NP/' % os.getcwd(),

        )

    home = os.path.expanduser("~")
    beastpath = '%s/BEASTv1.8.1/bin/beast' % home
    time = 130

    for folder in folders:
        
        xmlfile = '%sbeastinfile.xml' % folder
        command = beastpath + ' ' + xmlfile

        sbatchfname = '%srun.sbatch' % folder
        makesbatchfile = True
        if makesbatchfile:
            sbatchfile = CreateSbatchFile(sbatchfname,time,command)
        runsbatch = True
        if runsbatch:
            cdir = os.chdir(folder)
            sbatch = os.system('sbatch run.sbatch')
            


main()