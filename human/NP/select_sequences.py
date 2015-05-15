''' This script selects lineage for a set of DNA sequences for human NP

Written by Heather Machkovech

Functions defined in this module
------------------------------------

*``selectlineage`` : returns a list of sequences containing desired lineage based on year and subtype from a FASTA file
*``selectforlength`` : returns a list of sequences filtered by length

Input files
------------
*``sequences.fasta`` : unparsed NP DNA sequences

Output files
-------------
*``parsed.fasta`` : parsed NP sequences
'''

import listFASTA
import random
import os
from Bio import SeqIO

def selectlineage(seqs):
  """This function selects lineages to only include if year is known and subtype/year combo is:

        - H1N1 1918 to 1957

        - H2N2 1957 to 1968
        
        - H3N2 1968 to 2013 

      This function requires that the year be contained in the first 4 characters of the title line. 
      It also requires the subtype to be contained in the last 4 characters of title line.

      This function randomizes the sequence set.
  """
  sequences = listFASTA.listFASTA(seqs)

  seed = 1
  random.seed(seed)
  random.shuffle(sequences)
  
  basename = os.path.splitext(seqs)[0] # use basename as filename for new fasta
  print "Selecting lineage for %s" % basename
  list_sequences = []

  years_to_retain = { # keyed by subtype, values are all years to retain
        'H1N1':[year for year in range(1918, 1958)],
        'H2N2':[year for year in range(1957, 1969)],
        'H3N2':[year for year in range(1968, 2014)],
        }
  nretained = ntotal = 0
  for entries in sequences:
   # print entries
    ntotal += 1
    year = str(entries.description)[0:4]
    subtype = str(entries.description)[-4 : ]
    if year != "unkn" and year != "___c":
      year = int(year)
    else:
      continue
    if (subtype in years_to_retain) and (year in years_to_retain[subtype]):
      list_sequences.append(entries)
      nretained += 1
    else:
      print "not retaining %s" % entries.description
  print "Retained %d of %d sequences as being in the correct lineage" % (nretained, ntotal)
  return list_sequences

def selectforlength(seqs):
  """This function removes sequences that are not at least 1345 long"""

  print "Selecting sequences that are at least 1345 nt long" 
  nretained = ntotal = 0 
  lsequences = []
  for entries in seqs:
    ntotal += 1
    if len(entries.seq) >= 1345:
      nretained += 1
      print "%s" % entries.description, len(entries.seq)
      lsequences.append(entries)
    elif len(entries.seq) < 1345:
      print "not including %s" % entries.description
      print len(entries.seq)
    else:
      print "length not calculated for %s" % entries.description
  print 'Retained %d of %d sequences as being the correct length' % (nretained, ntotal)
 # print len(lsequences)
  return lsequences

def main():
  slineage = True
  selectlength = True

  if slineage:
    sequenceset = '%s/human/NP/sequences.fasta' % os.getcwd()
    sequences = selectlineage(sequenceset)
  if selectlength:
    sequences = selectforlength(sequences)
    SeqIO.write(sequences, "%s/human/NP/parsed.fasta"  % os.getcwd(), "fasta")

main()



