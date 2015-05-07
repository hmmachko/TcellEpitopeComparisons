def listFASTA(name):
	""" This function uses biopython to read a file *name* that is in fasta format and 
    returns a list of the sequence entries. 

    Each entry contains fields such as seq that contains the sequence and description that contains the 
    first line of the fasta file. See http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec55 
    """
	import os	

    
	#name = []
	basename = os.path.splitext(name)[0]
#	print basename
	records = '%s' % basename
	print "returning fasta list from ", records
	from Bio import SeqIO
	handle = open(name, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	print len(records)
	return records

	handle.close() 

