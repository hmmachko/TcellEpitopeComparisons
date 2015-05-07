import os

def makeepitopefile(infile,outfile):
    f = open(infile,'r')
    fx = open(outfile,'w')
    fx.write('site,epitopes\n')
    for lines in f:
        entry = lines.strip().split()
       # print entry
        
        for index in range(len(entry[0])):
           # print '%s %s' %(index, entry[0][index])
            fx.write('%s,%s\n' %(index+1,entry[0][index]))
    f.close()
    fx.close()

def main():
    inputepitopesites = '%s/human/HA_H3/h3-epitope-mask.tsv' % os.getcwd()
    outputepitopesites = '%s/human/HA_H3/antibodyepitopesbysite.csv' % os.getcwd()
    formatfile = True
    if formatfile:
        epitopefile = makeepitopefile(inputepitopesites,outputepitopesites)


main()

