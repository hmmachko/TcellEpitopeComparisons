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

    combinedRAxML = False
    combRAxML = '%s/scripts/raxmlhumanswinecombtree.py' % os.getcwd()
    if combinedRAxML:
        RAxML = subprocess.call(['python', combRAxML])

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
    avemutationpersite = '%s/scripts/runcalcavesub.sbatch' % os.getcwd()
    if calculateavemutationpersite:
        mutationspersite = subprocess.call(['sbatch',avemutationpersite])

    calcmutationantigentononantigen = False
    antigenictononantigeniccalc = '%s/scripts/calcavemutationantigenicnonantigenic.py' % os.getcwd()
    if calcmutationantigentononantigen:        
        antigenvsnonantigen = subprocess.call(['python',antigenictononantigeniccalc])

    plotantigentononantigen = True
    plotantigenictononantigenic = '%s/scripts/plotavemutationepitopetononepitope.py' % os.getcwd()
    if plotantigentononantigen:
        plotantigenvsnonantigen =  subprocess.call(['python',plotantigenictononantigenic])

    makenexus=False
    nexus = '%s/scripts/simplifytree.py' % os.getcwd()
    if makenexus:
        nexusfile = subprocess.call(['python',nexus])

    analyzednds = False
    dnds = '%s/scripts/dndsanalysis.py' % os.getcwd()
    if analyzednds:
        boxplot = subprocess.call(['python',dnds])
    calcf = False
    fcalculation = '%s/scripts/fanalysis.py' % os.getcwd()
    if calcf:
        f = subprocess.call(['python',fcalculation])

    plotf = True
    fplot = '%s/scripts/fplots.py' % os.getcwd()
    if plotf:
        f = subprocess.call(['python2',fplot])

    maxdivergence = False
    divergencecalc = '%s/scripts/maxsequencedivergence.py' % os.getcwd()
    if maxdivergence:
        mdiverg = subprocess.call(['python',divergencecalc])

    trunksubstitutionsummary = True
    trunksubsum = '%s/scripts/runavesubtimeidentity.sbatch' % os.getcwd()
    if trunksubstitutionsummary:
        tsum = subprocess.call(['sbatch',trunksubsum])

    #part2 = False
    #subsummary = '%s/scripts/trunksubstitutionidentityandtime.py' % os.getcwd()
    #if part2:
        #tsum = subprocess.call(['python',subsummary])


main()