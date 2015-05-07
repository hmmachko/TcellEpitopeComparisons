
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
  dnaseqmatch = re.compile('^[ACTGactg-]+$')
  sequences = listFASTA.listFASTA(seqs)
  for entries in sequences:
    ntotal += 1
    if dnaseqmatch.search(str(entries.seq)):
      list_sequences.append(entries)
      nretained += 1
    else:
      print "not including:", entries.description
      print entries.seq
  print 'Retained %d of %d sequences as not containing ambiguous characters ' % (nretained, ntotal)
  return list_sequences


def norepeats(seqs):
  """ This function gets rid of repeats in a list of sequences"""

  print "Selecting non-repeated sequences"

  nretained = ntotal = 0
  sequence_dict = {}
  for entry in seqs:
    entry.description = re.sub("[':)(;,]", "", entry.description)
    entry.description = re.sub(' ','',entry.description)
    entry.description = re.sub("'","",entry.description)
    entry.id = ''
    print entry.description
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
  namelist = []
  if host =='swine':
    subtype_refine = {'H3N2':{}}
  if host =='human':
    subtype_refine = {'H3N2':{}}
  sequences = []
  n = 3 #number of desired files per year/subtype
  nretained = ntotal = 0
  for entries in seqs:
    print entries.description
    ntotal += 1
    year = str(entries.description)[0:4]
    subt = 'H3N2'
    if year not in subtype_refine[subt]:
      #if len(entries.seq) >= (.9*1701):
      lengthprot = re.findall('[ACTG]',str(entries.seq))
      print 'checkingproteinlength'
      print len(lengthprot)
      if len(lengthprot) >= (.9*1701):

        sequences.append(entries)
        subtype_refine[subt][year] = 1
        nretained += 1
        namelist.append(entries.description)
    elif year in subtype_refine[subt]:
      if subtype_refine[subt][year] < n:
        lengthprot = re.findall('[ACTG]',str(entries.seq))
        print 'checkingproteinlength'
        print len(lengthprot)
        if len(lengthprot) >= (.9*1701):
        #if len(entries.seq) >= (.9*1701):
          nretained += 1
          sequences.append(entries)
          subtype_refine[subt][year] += 1
          namelist.append(entries.description)
      elif subtype_refine[subt][year] >= n:
        continue
  for entries in seqs:
    print entries.description
    
    year = str(entries.description)[0:4]
    subt = 'H3N2'
    if year not in subtype_refine[subt]:
      sequences.append(entries)
      subtype_refine[subt][year] = 1
      nretained += 1
    elif year in subtype_refine[subt]:
      if subtype_refine[subt][year] < n:
        if entries.description not in namelist:
        
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
  raxml_masterdir = '%sRAxML_output/' % basename

  raxml_dir = '%sRAxML_output/individual/' % basename
  treefile = 'prelimTree'
  raxml_tree = "%sRAxML_bestTree.%s" % (raxml_dir, treefile)
  if not os.path.isdir(raxml_masterdir):
    os.mkdir(raxml_masterdir)
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

def make_headers_sequences(sequences,host):
  header_seq = []
  for seq in sequences:
    print len(seq.seq)
    print len(seq.seq[:-2])
    if host =='swine':
      newseq = seq.seq[:-29]
    else: 
      newseq = seq.seq[:-3]
    #seq.description = re.sub("[':)(;,]", "", seq.description)
    #seq.description = re.sub(' ','',seq.description)
    h_s = (seq.description, newseq)
    header_seq.append(h_s)
  return header_seq


def Translateseqs(headers_sequences, readthrough_n=False, readthrough_stop=False, truncate_incomplete=False, translate_gaps=False):
    """Translates a set of nucleotide sequences to amino acid sequences.
    CALLING VARIABLES:
    `headers_sequences` : list of tuples `(header, seq)` as would be returned
    by `Read`.  The sequences should all specify valid coding nucleotide
    sequences.  
    
    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.
    `readthrough_n` : specifies that if any nucleotides
    in the sequence are equal to to an ambiguous nt code and cannot therefore
    be unambiguously translated into an amino acid, we simply translate through these 
    nucleotides by making the corresponding amino acid equal to "X".  By
    default, this option is `False`.  Note that even when this option is `False`,
    certain ambiguous nucleotides may still be translatable if they all lead to 
    the same amino acid.
    `readthrough_stop` : specifies that if we encounter any stop
    we simply translation them to 'X'.  By default,
    this option is `False`, meaning that we instead raise an error
    of an incomplete stop codon.
    `truncate_incomplete` : specifies that if the sequence
    length is not a multiple of three, we simply truncate off the one or two
    final nucleotides to make the length a multiple of three prior to translation.
    By default, this option is `False`, meaning that no such truncation is done.
    `translate_gaps` : specifies that a codon of '---' is translated
    to '-'. Codons with one '-' are also translated to gaps.
    RETURN VARIABLE:
    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.
    EXAMPLES:
    >>> Translate([('seq1', 'ATGTAA'), ('seq2', 'gggtgc')])
    [('seq1', 'M'), ('seq2', 'GC')]
    >>> Translate([('seq2', 'GGNTGC')])
    [('seq2', 'GC')]
    
    >>> Translate([('seq2', 'NGGTGC')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate codon NGG
    
    >>> Translate([('seq2', 'NGGTGC')], readthrough_n=True)
    [('seq2', 'XC')]
    >>> Translate([('seq2', 'TAATGC')])
    Traceback (most recent call last):
       ...
    ValueError: Premature stop codon
    >>> Translate([('seq2', 'TAATGC')], readthrough_stop=True)
    [('seq2', 'XC')]
    >>> Translate([('seq2', 'TGCA')])
    Traceback (most recent call last):
       ...
    ValueError: Sequence length is not a multiple of three
    >>> Translate([('seq2', 'TGCA')], truncate_incomplete=True)
    [('seq2', 'C')]
    >>> Translate([('seq2', 'TGC---')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate gap.
    >>> Translate([('seq2', 'TGC---')], translate_gaps=True)
    [('seq2', 'C-')]
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    assert isinstance(headers_sequences, list)
    translated_headers_sequences = []
    for (head, seq) in headers_sequences:
        seq = seq.upper()
        print len(seq)
        if len(seq) % 3:
            if truncate_incomplete:
                seq = seq[ : -(len(seq) % 3)]
            else:
                raise ValueError, "Sequence length is not a multiple of three"
        prot_length = len(seq) // 3
        prot = []
        for i in range(prot_length):
            codon = seq[3 * i : 3 * (i + 1)]
            try:
                aa = genetic_code[codon]
            except KeyError:
                if '-' in codon:
                    if translate_gaps:
                        aa = '-'
                    else:
                        raise ValueError("Cannot translate gap.")
                else:
                    # see if we have an ambiguous nucleotide codon that doesn't matter in the translation
                    possible_nt1 = AmbiguousNTCodes(codon[0])
                    possible_nt2 = AmbiguousNTCodes(codon[1])
                    possible_nt3 = AmbiguousNTCodes(codon[2])
                    possible_codons = []
                    for nt1 in possible_nt1:
                        for nt2 in possible_nt2:
                            for nt3 in possible_nt3:
                                possible_codons.append("%s%s%s" % (nt1, nt2, nt3))
                    try:
                        aa = genetic_code[possible_codons[0]]
                    except KeyError:
                        raise KeyError("Cannot translate codon %s in %s" % (codon, head))
                    for possible_codon in possible_codons:
                        if genetic_code[possible_codon] != aa:
                            if readthrough_n:
                                aa = 'X'
                            else:
                                raise ValueError("Cannot translate codon %s" % codon)
            if aa == 'STOP' and i == prot_length - 1:
                aa = ''
            elif aa == 'STOP':
                if readthrough_stop:
                    aa = 'X'
                else:
                    raise ValueError("Premature stop codon")
            prot.append(aa)
        translated_headers_sequences.append((head, ''.join(prot)))
    return translated_headers_sequences


def AmbiguousNTCodes(nt):
    """Returns all possible nucleotides corresponding to an ambiguous code.
    This method takes as input a single nucleotide character, which is
    assumed to represent a nucleotide as one of the accepted
    codes for an ambiguous character.  Returns a list giving
    all possible codes for which a nucleotide might stand.  Raises
    an exception if `nt` is not a valid nucleotide code.
    EXAMPLES:
    >>> AmbiguousNTCodes('N')
    ['A', 'T', 'G', 'C']
    >>> AmbiguousNTCodes('R')
    ['A', 'G']
    >>> AmbiguousNTCodes('A')
    ['A']
    >>> AmbiguousNTCodes('-')
    ['-']
    >>> AmbiguousNTCodes('F')
    Traceback (most recent call last):
       ...
    ValueError: Invalid nt code of "F"
    """
    if nt in ['A', 'T', 'G', 'C', '-']:
        return [nt]
    elif nt == 'R':
        return ['A', 'G']
    elif nt == 'Y':
        return ['T', 'C']
    elif nt == 'K':
        return ['G', 'T']
    elif nt == 'M':
        return ['A', 'C']
    elif nt == 'S':
        return ['G', 'C']
    elif nt == 'W':
        return ['A', 'T']
    elif nt == 'B':
        return ['C', 'G', 'T']
    elif nt == 'D':
        return ['A', 'G', 'T']
    elif nt == 'H':
        return ['A', 'C', 'T']
    elif nt == 'V':
        return ['A', 'C', 'G']
    elif nt == 'N':
        return ['A', 'T', 'G', 'C']
    else: 
        raise ValueError('Invalid nt code of "%s"' % nt)

def writefasta(sequences, outfile):
  f = open(outfile,'w')
  for seq in sequences:
    f.write('>%s\n'%seq[0])
    if seq[1][-1] == '-':
      f.write('%s\n'%seq[1][:-1])
    else:
      f.write('%s\n'%seq[1])
  f.close

def main():

  hosts = ('swine','human')


  for host in hosts:
    starting_sequencesets = (
   
    '%s/%s/HA_H3/select_sequences.py' % (os.getcwd(), host),     
 
    )

  
  
    for seqsets in starting_sequencesets:
      base = seqsets[:-19]
      print base

      select = True
      remove_ambiguous = True
      remove_anomalous = True
      remove_repeats = True
      threesequencespersubtypeperyear = True
      formatfortranslation = False
      translate_seq = False
      record_translated_seq = False
      alignment = False
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
        print "checking"
        print len(sequences)
        f = open("%scds_aligned.fasta" % base,'w')
        for seq in sequences:
          #f = open("%scds_aligned.fasta" % base,'w')
          seq.description = re.sub("[':)(;,]", "", seq.description)
          if host == 'swine':
            shortenedseq = seq.seq[:-32]
          if host =='human':
            shortenedseq = seq.seq[:-3]
          f.write('>%s\n' % seq.description)
          f.write('%s\n'%shortenedseq)
          print seq.description
        f.close()
        #SeqIO.write(sequences, "%scds_aligned.fasta" % base, "fasta")
      if formatfortranslation:
        heads_seqs = make_headers_sequences(sequences,host)
      if translate_seq:
        translatedseq = Translateseqs(heads_seqs,translate_gaps=True)
        for seqs in translatedseq:
          print seqs
      if record_translated_seq:
        proteinoutfile = '%sprot_aligned.fasta' % base
        translatedrecord = writefasta(translatedseq, proteinoutfile)
        #seqs_for_trans = ('%sparsed.fasta' % base,)
        #for seqset in seqs_for_trans:
          #translate_seq = translate(seqset)

      if alignment:
        reference = ('%sreference_sequence_translated.fasta' % base)
        proteins = ('%sparsed_translated.fasta' % base)
        alignment = align(reference, proteins)
      if buildtree:
        sequencesets = '%scds_aligned.fasta' % base
        tree = raxmltree(sequencesets)



main()