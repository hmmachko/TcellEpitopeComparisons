import os

def readsites(sitefile,intype):
    site_dict = {} 

    site_file = open(sitefile, 'r')
    with site_file as f:

        #if sitefile == ('%s/Human/HA_H1/caton_epitopes_by_site2.csv' % os.getcwd()):
        next(f)
        for lines in site_file:     
            entry = lines.strip().split(",")
            site_dict[int(entry[0])] = intype(entry[1])
    site_file.close()
    #print sitefile
    #print site_dict
    return site_dict

# calc ave mutation for epitopes vs non epitopes = number of mutations for ep/number of epitopes
def startoutput(output):
    with open(output, "w") as f:
        f.write("ratio of epitope to nonepitope mutations\n")
        f.close()
def calcavemutation(sitedict, mutationdict,outlabel,output,protein): 
    count_epitope = 0
    epitope_mutations = 0
    count_nonepitope = 0
    nonepitope_mutations = 0  
            
    for aa in sitedict:
       
        if sitedict[aa] != 0:
            count_epitope +=1
            #print mutationdict[aa]
            epitope_mutations += mutationdict[aa]
            #print aa, mutationdict[aa], count_epitope
        elif sitedict[aa] ==0:
            count_nonepitope +=1
            #print mutationdict[aa]

            nonepitope_mutations += mutationdict[aa]
        #else:
            #print "non 1/0 number"
    ave_epitope_mutation = float(epitope_mutations)/count_epitope
    ave_nonepitope_mutation = float(nonepitope_mutations)/count_nonepitope
    ratio_epitope_to_nonepitope = float(ave_epitope_mutation)/ave_nonepitope_mutation
    print "%s,%s" % (outlabel, ratio_epitope_to_nonepitope)

    with open(output, "a") as f:
        #f.write("%s\n"% protein)
        f.write("%s %s,%s,%s,%s,%s,%s,%s,%s\n" % (protein,outlabel, ratio_epitope_to_nonepitope,epitope_mutations,count_epitope,ave_epitope_mutation, nonepitope_mutations,count_nonepitope,ave_nonepitope_mutation))
        f.close()
   



def main():
    epitopetypes = ('antibody',)
    
    for epitopetype in epitopetypes:
        if epitopetype == 'antibody':
            folder = 'antibody'
        else:
            folder = epitopetype
        outfile_name = ('%s/plots/%s/avemutationratio_epitope_nonepitope.csv'% (os.getcwd(), folder))
        startoutfile = True
        if startoutfile:
            outf = startoutput(outfile_name)
        
        if epitopetype == 'antibody':         
            proteins = ('HA_H3',)

        if epitopetype == 'cd8':
            proteins = ('M1','NP')

        for protein in proteins:
            with open(outfile_name, "a") as f:
                f.write('%s, average_mutation_epitope_to_nonepitope,number_epitopemutations,number_epitopesites, numerator,number_nonepitopemutations,number_nonepitopesites,denominator,\n'%protein)
                f.close()
            print protein
            if epitopetype == 'cd8':
                epitopesbysite = '%s/human/%s/%scombinedepitopesbysite.csv' % (os.getcwd(), protein,epitopetype)
            if epitopetype == 'antibody':
                epitopesbysite= ('%s/human/HA_H3/antibodyepitopesbysite.csv' % os.getcwd())
            avemutationpersite = (
                        '%s/human/%s/trunkavemutationpersite.csv'% (os.getcwd(),protein),
                        '%s/swine/%s/trunkavemutationpersite.csv'% (os.getcwd(),protein),
                        '%s/human/%s/treeavemutationpersite.csv'% (os.getcwd(),protein),
                        '%s/swine/%s/treeavemutationpersite.csv'% (os.getcwd(),protein),
                        )
            outfile_labels = ['human trunk', 'swine trunk', 'human tree', 'swine tree']
            file_and_label = list(zip(avemutationpersite,outfile_labels))
      
            epitopestodict = True
            if epitopestodict:
                epitopesites = readsites(epitopesbysite,int)
 
            for files in file_and_label:
                avemutationspersitetodict = True
                if avemutationspersitetodict:
                    mutationstodict = readsites(files[0],float)
                avemutation = True
                if avemutation:
                    mutation = calcavemutation(epitopesites, mutationstodict,files[1],outfile_name,protein)





main()




