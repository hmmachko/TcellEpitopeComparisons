'''This script calculates the average substitution rate for epitope to non-epitope regions. This is 
done for human and swine NP and M1. This is done for the entire tree and from the trunk only. The result 
for all proteins and trees is written to a csv file.

Infiles
---------
*``combinedepitopesbysite.csv`` : file with number of epitopes per site
*``trunkavemutationpersite.csv`` : file with average number of mutations per site for the trunk
*``treeavemutationpersite.csv`` : file with average number of mutations per site for the tree

Outfiles
-----------
*``avemutationratio_epitope_nonepitope.csv`` : summary of mutation rate of epitope to non-epitope regions, in plots/ directory

'''

import os

def readsites(sitefile,intype):
    '''This funciton reads a file *sitefile* that has a titleline, and subsequent lines with amino acid
    number at the first entry and another integer or float second entry (given by intype -float or int). 
    This function returns a dictionary with amino acid site as key and second entry as value. 
    '''
    site_dict = {} 

    site_file = open(sitefile, 'r')
    with site_file as f:

        next(f)
        for lines in site_file:     
            entry = lines.strip().split(",")
            site_dict[int(entry[0])] = intype(entry[1])
    site_file.close()
    return site_dict

def startoutput(output):
    '''This function initiates the *output* outfile for the average mutation rate in epitopes to non-epitopes.'''
    with open(output, "w") as f:
        f.write("ratio of epitope to nonepitope mutations\n")
        f.close()

def calcavemutation(sitedict, mutationdict,outlabel,output,protein): 
    '''This function calculates the average mutation rate in epitopes to non-epitopes. 
    *sitedict* is the dictionary that contains the amino acid site as key and number of epitopes as value.
    *mutationdict* is the dictionary that contains the amino acid site as key and average number of mutations as value.
    *outlabel* string indicating tree type of analysis, ex human trunk
    *output* is the outfile containing the average mutation rate 
    *protein* protein
    '''
    count_epitope = 0
    epitope_mutations = 0
    count_nonepitope = 0
    nonepitope_mutations = 0  
            
    for aa in sitedict:
       
        if sitedict[aa] != 0:
            count_epitope +=1
            epitope_mutations += mutationdict[aa]
        elif sitedict[aa] ==0:
            count_nonepitope +=1
            nonepitope_mutations += mutationdict[aa]
    ave_epitope_mutation = float(epitope_mutations)/count_epitope
    ave_nonepitope_mutation = float(nonepitope_mutations)/count_nonepitope
    ratio_epitope_to_nonepitope = float(ave_epitope_mutation)/ave_nonepitope_mutation
    print "%s,%s" % (outlabel, ratio_epitope_to_nonepitope)

    with open(output, "a") as f:
        #f.write("%s\n"% protein)
        f.write("%s %s,%s,%s,%s,%s,%s,%s,%s\n" % (protein,outlabel, ratio_epitope_to_nonepitope,epitope_mutations,count_epitope,ave_epitope_mutation, nonepitope_mutations,count_nonepitope,ave_nonepitope_mutation))
        f.close()
   
def main():
    epitopetypes = ('cd8',)
    
    for epitopetype in epitopetypes:
        
        folder = epitopetype
        outfile_name = ('%s/plots/%s/avemutationratio_epitope_nonepitope.csv'% (os.getcwd(), folder))
        startoutfile = True
        if startoutfile:
            outf = startoutput(outfile_name)
        proteins = ('M1','NP')

        for protein in proteins:
            with open(outfile_name, "a") as f:
                f.write('%s, average_mutation_epitope_to_nonepitope,number_epitopemutations,number_epitopesites, numerator,number_nonepitopemutations,number_nonepitopesites,denominator,\n'%protein)
                f.close()
            print protein
            
            epitopesbysite = '%s/human/%s/%scombinedepitopesbysite.csv' % (os.getcwd(), protein,epitopetype)
           
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




