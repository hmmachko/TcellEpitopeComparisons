'''This script submits for sbatch the command to build BEAST trees for influenza NP, PA, M1, and H1 HA. 

Functions
-----------
*``CreateSbatchFile`` : creates sbatch file

Input files
-------------
*``beastinfile.xml`` : input file to run BEAST

Output files
--------------
*``run.sbatch`` : sbatch file to run BEAST
*``prot_aligned.trees`` : BEAST output file containing thinned trees 
*``prot_aligned.log` : BEAST log output

Written by Heather Machkovech Feb 2015

'''

import os
import sys

def CreateSbatchFile(fname,time,command):
    f = open(fname,'w')
    f.write('#!/bin/sh\n')
    f.write('#SBATCH\n')
    f.write('#PBS -l walltime=%s:00:00\n' % time)
    f.write('%s'%command)
    f.close()


def main():

    folders = (
        #'%s/human/M1/' % os.getcwd(),
        #'%s/human/NP/' % os.getcwd(),
        #'%s/swine/M1/' % os.getcwd(),
       # '%s/swine/NP/' % os.getcwd(),
        '%s/human/HA_H3/' % os.getcwd(),
        '%s/swine/HA_H3/' % os.getcwd(),
        )

    home = os.path.expanduser("~")
    beastpath = '%s/BEASTv1.8.1/bin/beast' % home
    time = 100

    for folder in folders:
        xmlfile = '%sseparatebeastinfile.xml' % folder
        sbatchfname = '%sseparatetreerun.sbatch' % folder
        command = beastpath + ' ' + xmlfile

        makesbatchfile = True
        if makesbatchfile:
            sbatchfile = CreateSbatchFile(sbatchfname,time,command)
        runsbatch = True
        if runsbatch:
            cdir = os.chdir(folder)
            sbatch = os.system('sbatch separatetreerun.sbatch')
            


main()