
import subprocess
import os
import re
import listFASTA
from Bio import SeqIO
import random
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import glob


def noambiguous(seqs):
  """This function removes sequences that contain ambiguous sequence characters. """

  print "Selecting sequences without ambiguous characters"
  nretained = ntotal = 0 
  list_sequences = []
  dnaseqmatch = re.compile('^[ACTGactg]+$')
  sequences = listFASTA.listFASTA(seqs)
  for entries in sequences:
    ntotal += 1
    if dnaseqmatch.search(str(entries.seq)):
      list_sequences.append(entries)
      nretained += 1
    else:
      print "not including:", entries.description
      #print entries.seq
  print 'Retained %d of %d sequences as not containing ambiguous characters ' % (nretained, ntotal)
  return list_sequences


def norepeats(seqs):
  """ This function gets rid of repeats in a list of sequences"""

  print "Selecting non-repeated sequences"

  nretained = ntotal = 0
  sequence_dict = {}
  for entry in seqs:
    ntotal += 1
    stringentry = str(entry.seq)
    if stringentry not in sequence_dict:
      sequence_dict[stringentry] = entry
      nretained += 1
  unique_entries = sequence_dict.values() 
  print 'Retained %d of %d sequences as unique sequences ' % (nretained, ntotal)  
  return unique_entries
 

def noanomalies(seqs):
  """This function removes previously identified anomalous sequences that are 
  found in JVI_82_8947_Anomalies.txt, JDB_Anomalies.txt, and HMM_Anomalies.txt."""

  print "Removing anomalous sequences"
  sequences = []
  anomalies = frozenset([line.strip() for line in open('%s/inputfiles/anomalousseq/JVI_82_8947_Anomalies.txt'%os.getcwd()).readlines() + open('%s/inputfiles/anomalousseq/HMM_Anomalies.txt'%os.getcwd()).readlines() + open('%s/inputfiles/anomalousseq/JDB_Anomalies.txt'%os.getcwd()).readlines()]) # strain names of anomalous sequences 
  nretained = ntotal = 0
  for entries in seqs:
    ntotal += 1
    retain_entry = True
    for anomaly in anomalies:
      if anomaly in entries.description:
        retain_entry = False
        print "anomaly found %s:" % entries.description
        break
    if retain_entry:
      sequences.append(entries)
      nretained += 1
  print 'Retained %d of %d sequences as not containing anomalous sequences' % (nretained, ntotal)
  return sequences


def clean(seqs,host):
  """This function keeps 3 sequences/year per subtype. It requires that the year be contained 
  in the first 4 characters of the title line."""
  
  subtype_refine = {}
  if host =='swine':
    subtype_refine = {'H1N1':{}, 'H1N2':{}, 'H3N2':{}}
  if host =='human':
    subtype_refine = {'H1N1':{}, 'H2N2':{}, 'H3N2':{}}
  sequences = []
  n = 3 #number of desired files per year/subtype
  nretained = ntotal = 0
  for entries in seqs:
    ntotal += 1
    year = str(entries.description)[0:4]
    subt = str(entries.description)[-4 : ]
    if year not in subtype_refine[subt]:
      sequences.append(entries)
      subtype_refine[subt][year] = 1
      nretained += 1
    elif year in subtype_refine[subt]:
      if subtype_refine[subt][year] < n:
        nretained += 1
        sequences.append(entries)
        subtype_refine[subt][year] += 1
      elif subtype_refine[subt][year] >= n:
        continue

  print subtype_refine
  print 'Retained %d of %d sequences so that there up to 3 sequences per subtype per year' % (nretained, ntotal)
  return sequences
  
def make_protein_record(nuc_record):
  """This function returns a new SeqRecord with the translated sequence."""
  
  return SeqRecord(seq = nuc_record.seq.translate(cds=False), \
                    id = nuc_record.id, \
                    description = nuc_record.description)

def translate(sequences): 
  """This function translates sequences and writes output to a FASTA file."""
  
  basename = os.path.splitext(sequences)[0]
  print "Now translating sequences"
  proteins = (make_protein_record(nuc_rec) for nuc_rec in \
            SeqIO.parse(sequences, "fasta"))
  SeqIO.write(proteins, "%s_translated.fasta" % basename, "fasta")

def align(reference, proteins):
  """This function performs pairwise alignments of sequences in a fasta file and returns the sequences as both a protein and
  cDNA fasta file with the gaps stripped relative to reference. It also removes characters from the fasta defline not allowed in RAxML.
  It also subtracts 24 yrs from H1N1 from 1977 on. 

  reference = reference protein for pairwise alignments
  proteins = fasta sequences to be aligned

  """

  print "Now aligning sequences"
  basename = proteins[:-23]
 
  needle_cline = NeedleCommandline(asequence=reference, bsequence=proteins, gapopen=10, gapextend=0.5, outfile="%saligned.fasta" % basename, aformat='fasta')
  stdout, stderr = needle_cline()

  # the output alignment.fasta contains: reference, protein 1, reference, protein 2 etc, so now we want 
  # to take only protein 1, protein 2 etc, and remove the gaps relative to the reference

  align = listFASTA.listFASTA('%saligned.fasta' % basename)
  prots = listFASTA.listFASTA('%sparsed_translated.fasta' % basename)
  DNAseqs = listFASTA.listFASTA('%sparsed.fasta'% basename)
  refprot = listFASTA.listFASTA('%sreference_sequence_translated.fasta'% basename)

  refprotein = []
  for refpro in refprot:
    refprotein.append(refpro.seq)

# make list with desired defline >{year}/{month}/{day}{accession}{strain}{segname}{serotype} and remove ':)(;, (not allowed in RAxML) 
  heads = []
  seqs = []
  print "removing forbidden characters,adjusting year for H1N1 from 1977 on, and removing spaces in defline"
  for i in DNAseqs:
    seqs.append(i.seq)
    #print "description", i.description
    #print "id", i.id
    i.id = ''
    forbidden = re.sub("[':)(;,]", "", i.description)
    #print forbidden
    i.description = forbidden
    fspace = re.sub(' ', '', i.description)
    i.description = fspace
    year = int(i.description[0:4])
    subtype = str(i.description[-4 : ])
    #if (year >= 1977 and subtype == 'H1N1'):
      #adjustedyear = year - 24
      #base = str(i.description[4:])
      #adjustyear = re.sub(i.description[0:4], "%s" % newyear, i.description)
      #i.description = str(adjustedyear) + base
      #heads.append(i.description)
    #else:
    heads.append(i.description)

  assert len(align) == 2 * len(prots) == 2 * len(seqs)
  prot_alignment = []
  cds_alignment = []
  for i in range(len(prots)):
    prot = prots[i]
    head = heads[i]
    seq = seqs[i]
    
    (refa, prota) = (align[2 * i], align[2 * i + 1])

    assert len(refa) == len(prota)
    iref = iprot = 0
    alignedprot = []
    alignedcds = []
    for (aa_ref, aa_prot) in zip(refa.seq, prota.seq):
      assert (aa_ref == '-' or aa_ref == refpro[iref])
      assert (aa_prot == '-' or aa_prot == prot[iprot])
      if aa_ref == '-' and aa_prot != '-':
        iprot += 1
      elif aa_prot == '-' and aa_ref != '-':
        alignedprot.append(aa_prot)
        alignedcds.append('---')
        iref += 1
      elif aa_ref != '-' and aa_prot != '-':
        alignedprot.append(aa_prot)
        alignedcds.append(str(seq[3 * iprot : 3 * iprot + 3]))
        iref += 1
        iprot += 1
      else:
        raise ValueError("Both prots in alignment have gap")
    
    alignedprot = ''.join(alignedprot)
    alignedcds = ''.join(alignedcds)
    #print alignedprot
    #print alignedcds
    assert len(alignedprot) == len(refpro)
    protrecord = SeqRecord(Seq(alignedprot), id = head, description = '')
    prot_alignment.append(protrecord)
    cdsrecord = SeqRecord(Seq(alignedcds), id = head, description = '')
    cds_alignment.append(cdsrecord)
  assert len(prot_alignment) == len(cds_alignment)
  SeqIO.write(prot_alignment, "%sprot_aligned.fasta" % basename, "fasta")
  SeqIO.write(cds_alignment, "%scds_aligned.fasta" % basename, "fasta")

def raxmltree(sequencesets):
  """ This function generates a quick RAxML tree to visualize anomalous sequences that do not fit the molecular clock.
  sequencesets = fasta file containing sequences for RAxML tree

  """
  use_existing_output = False
  basename = sequencesets[:-17]
  raxml_dir = '%sRAxML_output/individual/' % basename
  treefile = 'prelimTree'
  raxml_tree = "%sRAxML_bestTree.%s" % (raxml_dir, treefile)
  if not os.path.isdir(raxml_dir):
    os.mkdir(raxml_dir)
  print "Using RAxML to build a quick phylogenetic tree for visual analysis of potentially anomalous sequences..." 
  print "\nFor sequence set %s" % sequencesets
  if use_existing_output and os.path.isfile(raxml_tree):
    print "The existing tree of %s will be used, and no new tree will be created." % raxml_tree
  else:
    oldfiles = glob.glob("%s/*%s*" % (raxml_dir, treefile))
    for fname in oldfiles:
      os.remove(fname)
    os.system('%s -w %s -n %s -p 1 -m GTRCAT -s %s --no-bfgs' % ('raxmlHPC-SSE3', raxml_dir, treefile, sequencesets))
  if not os.path.isfile(raxml_tree):
    raise ValueError("Failed to generated expected output tree file %s" % raxml_tree)
  print "RAxML tree complete"


def main():

  hosts = ('swine',)#'human',


  for host in hosts:
    starting_sequencesets = (
    '%s/%s/NS1/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/NS2/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/M1/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/M2/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/NP/select_sequences.py' % (os.getcwd(), host),  
    '%s/%s/HA_H1/select_sequences.py' % (os.getcwd(), host),     
    '%s/%s/NA_N1/select_sequences.py' % (os.getcwd(), host), 
    '%s/%s/PA/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/PB1_P1/select_sequences.py' % (os.getcwd(), host),
    '%s/%s/PB2/select_sequences.py' % (os.getcwd(), host),
    )
  
  
    for seqsets in starting_sequencesets:
      base = seqsets[:-19]
      print base

      select = True
      remove_ambiguous = True
      remove_anomalous = True
      remove_repeats = True
      threesequencespersubtypeperyear = True
      translate_seq = True
      alignment = True
      buildtree = True
    
      if select:
  	   sequences = subprocess.call(['python', seqsets])
  
      if remove_ambiguous:
        sequences = noambiguous('%sparsed.fasta' % base)
      if remove_anomalous:
        sequences = noanomalies(sequences)
      if remove_repeats:
        sequences = norepeats(sequences)
      if threesequencespersubtypeperyear:
        sequences = clean(sequences,host)
        print 'checking'
        for sequence in sequences:
          print sequence.description
        SeqIO.write(sequences, "%sparsed.fasta" % base, "fasta")
      if translate_seq:
        seqs_for_trans = ('%sparsed.fasta' % base, "%sreference_sequence.fasta" % base)
        for seqset in seqs_for_trans:
          translate_seq = translate(seqset)

      if alignment:
        reference = ('%sreference_sequence_translated.fasta' % base)
        proteins = ('%sparsed_translated.fasta' % base)
        alignment = align(reference, proteins)
      if buildtree:
        sequencesets = '%scds_aligned.fasta' % base
        tree = raxmltree(sequencesets)



main()