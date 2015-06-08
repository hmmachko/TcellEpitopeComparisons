import os
import listFASTA

def calcdivergence(inputfile,summary):
    sequencelist = listFASTA.listFASTA(inputfile)
    #f = open(formattedfile,'w')
    maxdiff = 0
    maxindexa = ''
    maxindexb = ''
    numberseq = len(sequencelist)
    proteinlength = len(sequencelist[0].seq)
    for index,seqs in enumerate(sequencelist):
        #print index
        bindex= index+1    
        for i,seqb in enumerate(sequencelist[bindex:]):
            numberdif = 0
            for sitenumber,site in enumerate(seqs.seq):
                #print sitenumber,site
                #print seqb
                if site != seqb.seq[sitenumber]:
                    numberdif +=1
                if sitenumber == proteinlength - 1:
                    if numberdif > maxdiff:
                        maxdiff = numberdif
                        maxindexa = seqs.description
                        maxindexb = seqb.description
    print inputfile
    print maxdiff
    diver = float(maxdiff)/proteinlength
    percentdiver = diver*100
    print 'percent divergence: %s' % percentdiver
    print maxindexa
    print maxindexb
    f = open(summary, 'a')
    f.write('%s\n' % inputfile)
    f.write('percent divergence: %s\n' % percentdiver)
    f.write('%s\n' % maxindexa)
    f.write('%s\n' % maxindexb)

    #f.close()

def main():
    fastafiles = [
        '%s/human/M1/prot_aligned.fasta' % os.getcwd(),
        '%s/human/NP/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/M1/prot_aligned.fasta' % os.getcwd(),
        '%s/swine/NP/prot_aligned.fasta' % os.getcwd(),]
    csvfile = [
        '%s/human/M1/prot_aligned_format.fasta' % os.getcwd(),
        '%s/human/NP/prot_aligned_format.fasta' % os.getcwd(),
        '%s/swine/M1/prot_aligned_format.fasta' % os.getcwd(),
        '%s/swine/NP/prot_aligned_format.fasta' % os.getcwd(),]
    summaryfile = '%s/plots/summarydivergence.csv' % os.getcwd()
    f = open(summaryfile, 'w')
    f.close()
    for seqfile in fastafiles:
        calcseqdivergence =True
        if calcseqdivergence:
            divergence = calcdivergence(seqfile,summaryfile)

main()