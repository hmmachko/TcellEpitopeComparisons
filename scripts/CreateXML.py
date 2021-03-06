'''
This script constructs BEAST input XML files in order to build phylogenic trees using BEAST for influenza NP and M1.

Functions
-----------
*``CreateHead`` : function that creates top portion of XML file that contains taxa names and sequences
*``JoinHeadTail`` : function that concatenates top and bottom portions of XML file

Input files
-------------
*``combined_prot_aligned.fasta`` : 1 file for each protein that contains aligned human and swine influenza protein sequences used for building trees
*``tail.xml`` : bottom portion of XML that contains parameters for BEAST

Output files
-------------
*``head.xml`` : 1 file generated for each protein containing top portion of XML file
*``beastinfile.xml`` : 1 complete XML file generated for each protein 

Written by Heather Machkovech Feb 2015
'''
import listFASTA
import os

def CreateHead(seqs, head):
    '''This function takes aligned protein sequences *seqs* and creates the top portion of the XML file saved as *head*.
    '''
    f = open(head, 'w') 
    f.write('<?xml version="1.0" standalone="yes"?>\n\n\
<!-- Generated by BEAUTi v1.8.1                                              -->\n\
<!--       by Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard         -->\n\
<!--       Department of Computer Science, University of Auckland and        -->\n\
<!--       Institute of Evolutionary Biology, University of Edinburgh        -->\n\
<!--       David Geffen School of Medicine, University of California, Los Angeles-->\n\
<!--       http://beast.bio.ed.ac.uk/                                        -->\n\
<beast>\n\n\
	<!-- The list of taxa to be analysed (can also include dates/ages).          -->\n')

    sequences = listFASTA.listFASTA(seqs)
    f.write('	<!-- ntax=%s      -->\n\
	<taxa id="taxa">\n' % len(sequences)) 
    for seqs in sequences:
        f.write('		<taxon id="%s">\n' % seqs.description)
        f.write('			<date value="%s.0" direction="forwards" units="years"/>\n' % seqs.description[0:4])
        f.write('		</taxon>\n')



	f.write('	</taxa>\n\n\
	<!-- The sequence alignment (each sequence refers to a taxon above).         -->\n\
	<!-- ntax=%s nchar=%s                                                      -->\n\
	<alignment id="alignment" dataType="amino acid">\n' % (len(sequences), len(seqs.seq)))

    print len(seqs.seq)
    for seqs in sequences:
		f.write('		<sequence>\n\
			<taxon idref="%s"/>\n\
			%s\n\
		</sequence>\n' % (seqs.description, seqs.seq))
    f.write('	</alignment>\n\n') 
    f.close()

def JoinHeadTail(head, tail, joined):
    '''This function concatenates two files (*head* to *tail*) and saves new file as *joined*.
    '''
    filenames = [head, tail]
    with open(joined, 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				outfile.write(infile.read())

def main():
    sequencesets = (
        '%s/human/M1/combined_prot_aligned.fasta' % os.getcwd(),
        '%s/human/NP/combined_prot_aligned.fasta' % os.getcwd(),  
       ) 

    for sequence in sequencesets:
        base = sequence[:-27]
        head = '%shead.xml' % base
        tail = '%s/inputfiles/combined_tail.xml' % os.getcwd()
        xmlinfile = '%sbeastinfile.xml' % base
    
        makehead = True
        if makehead:
            dochead = CreateHead(sequence, head)
        xmljoin = True
        if xmljoin:
            xml = os.system('cat %s %s > %s' % (head,tail,xmlinfile))


main()
	