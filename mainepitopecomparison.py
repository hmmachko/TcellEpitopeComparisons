import subprocess
import os

def main():

    SelectSeqs = False
    seqparser = '%s/scripts/parse_sequences.py' % os.getcwd()
    if SelectSeqs:
        parsedsequences = subprocess.call(['python', seqparser])

    combinedtree = False
    raxmltree = '%s/scripts/raxmlhumanswinecombtree.py' % os.getcwd()
    if combinedtree:
        combtree = subprocess.call(['python', raxmltree])

    sbatchepitopes = False
    epitoperunfile = '%s/scripts/sbatchepitopes.py' % os.getcwd()
    if sbatchepitopes:
        beast = subprocess.call(['python', epitoperunfile])

    plotepitopesperprotein = False
    epperprotein = '%s/scripts/plotepitopesperprotein.py' % os.getcwd()
    if plotepitopesperprotein:
        epitopesperprotein = subprocess.call(['python', epperprotein])

    makecombinedXMLfile=False
    xmlinputfile1 = '%s/scripts/CreateXML.py' % os.getcwd()
    if makecombinedXMLfile:
        xmlfile = subprocess.call(['python', xmlinputfile1])

    sbatchcombinedbeast = False
    beastrunfile1 = '%s/scripts/sbatchphylogenies.py' % os.getcwd()
    if sbatchcombinedbeast:
        beast = subprocess.call(['python', beastrunfile1])

    makeseparateXMLfile=False
    xmlinputfile2 = '%s/scripts/CreateSeparateXML.py' % os.getcwd()
    if makeseparateXMLfile:
        xmlfile = subprocess.call(['python', xmlinputfile2])

    sbatchseparatebeast = False
    beastrunfile2 = '%s/scripts/sbatchseparatephylogenies.py' % os.getcwd()
    if sbatchseparatebeast:
        beast = subprocess.call(['python', beastrunfile2])

    makemaxcladecredibilitytree = False
    maketree = '%s/scripts/MakeMaxCladeCred.py' % os.getcwd()
    if makemaxcladecredibilitytree:
        tree = subprocess.call(['python',maketree])

    calculateavemutationpersite = False
    avemutationpersite = '%s/scripts/calcaveragemutationpersite.py' % os.getcwd()
    if calculateavemutationpersite:
        mutationspersite = subprocess.call(['python',avemutationpersite])

    calcmutationantigentononantigen = False
    antigenictononantigeniccalc = '%s/scripts/calcavemutationantigenicnonantigenic.py' % os.getcwd()
    if calcmutationantigentononantigen:        
        antigenvsnonantigen = subprocess.call(['python',antigenictononantigeniccalc])

    plotantigentononantigen = False
    plotantigenictononantigenic = '%s/scripts/plotavemutationepitopetononepitope.py' % os.getcwd()
    if plotantigentononantigen:
        plotantigenvsnonantigen =  subprocess.call(['python',plotantigenictononantigenic])

    calcf = False
    fcalculation = '%s/scripts/fanalysis.py' % os.getcwd()
    if calcf:
        f = subprocess.call(['python',fcalculation])

    plotf = True
    fplot = '%s/scripts/fplots.py' % os.getcwd()
    if plotf:
        f = subprocess.call(['python2',fplot])

##


    calcavebranchmutation = False
    avebranchmut = '%s/scripts/calcbranchavemutation.py' % os.getcwd()
    if calcavebranchmutation:
        f = subprocess.call(['python',avebranchmut])

    calcfbranch = False
    fcalculationbranch = '%s/scripts/fanalysisbranch.py' % os.getcwd()
    if calcfbranch:
        f = subprocess.call(['python',fcalculationbranch])

    plotfbranch = False
    fplotbranch = '%s/scripts/fplotsbranch.py' % os.getcwd()
    if plotfbranch:
        f = subprocess.call(['python2',fplotbranch])


    dndsFEL = False
    FELboxplot = '%s/scripts/boxplotdnds.py' % os.getcwd()
    if dndsFEL:
        boxplot = subprocess.call(['python',FELboxplot])

    dndsFEL2 = False
    FELboxplot2 = '%s/scripts/FEL2.py' % os.getcwd()
    if dndsFEL2:
        print 'usertree'
        boxplot = subprocess.call(['python',FELboxplot2])

    analyzeFUBAR = False
    FUBARplots = '%s/scripts/FUBARanalysis.py' % os.getcwd()
    if analyzeFUBAR:
        boxplot = subprocess.call(['python2',FUBARplots])

    analyzeFUBAR2 = False
    FUBARplots2 = '%s/scripts/FUBARanalysis2.py' % os.getcwd()
    if analyzeFUBAR2:
        boxplot = subprocess.call(['python',FUBARplots2])
    makenexus = False
    nexus = '%s/scripts/simplifytree.py' % os.getcwd()
    if makenexus:
        nexusfile = subprocess.call(['python',nexus])

    selectHA = False
    
    selectHAseqs = '%s/scripts/parse_sequences_HA.py' % os.getcwd()
    if selectHA:
        parsedsequences = subprocess.call(['python', selectHAseqs])

    formatantibodyepitopes = False
    formatepitopes = '%s/scripts/remakeepitopesfile.py' % os.getcwd()
    if formatantibodyepitopes:
        epitopesfile = subprocess.call(['python', formatepitopes])









    # general steps:
        # make separate trees from combined trees
        # calc average mutations per site for whole tree and trunk
        # plot normalized epitope to non epitope mutation rate
        # calc 'f'


main()