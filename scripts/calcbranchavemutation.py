import os

def readfiletodictionary(filename, makelist = False):
    '''this function takes a ``filename`` where the file is in the following format:
    titleline
    1, 0
    2, 2
    3, 5
    and returns a dictionary with first entry as key and second entry as value
    There is an option to return list containing first column by setting makelist=True when function called
        '''

    fx = open(filename, 'r')
    sites = []
    filedict = {}
    with fx as f:
        next(f)
        for lines in fx:     
            entry = lines.strip().split(",")
            if makelist:           
                sites.append(int(entry[0]))
            filedict[int(entry[0])] = float(entry[1])
    fx.close()
    if makelist:
        return sites, filedict
    else:
        return filedict

def calcbranchavemutation(sitedicttree,sitedicttrunk,outfile):

    print 'starting'

    sitedictbranch = {}

    for aminoacid in sitedicttree:
        
        branchsubrate = sitedicttree[aminoacid] - sitedicttrunk[aminoacid]
        print '%s %s %s' % (sitedicttree[aminoacid],sitedicttrunk[aminoacid],branchsubrate)
        sitedictbranch[aminoacid] = branchsubrate

    f = open(outfile,'w')
    f.write('site,avemutation\n')
    for aminoacid in sorted(sitedictbranch):
        f.write('%s,%s\n' %(aminoacid, sitedictbranch[aminoacid]))

    f.close()

    return sitedictbranch

def main():
    proteins = ('M1','NP')
    hosts = ('swine','human')
    for protein in proteins:
        for host in hosts:
        
            treeavemutation = '%s/%s/%s/treeavemutationpersite.csv' % (os.getcwd(), host, protein)
            trunkavemutation = '%s/%s/%s/trunkavemutationpersite.csv' % (os.getcwd(), host,protein)
            branchavemutation = '%s/%s/%s/branchavemutationpersite.csv' % (os.getcwd(), host,protein)
            opentree = True
            if opentree:
                treemut = readfiletodictionary(treeavemutation)
                print treemut
            opentrunk = True
            if opentrunk:
                trunkmut = readfiletodictionary(trunkavemutation)
                print trunkmut
            getbranchmut = True
            if getbranchmut:
                branchmut = calcbranchavemutation(treemut,trunkmut,branchavemutation)

main()






