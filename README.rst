============================================
Analysis of CD8 T-cell epitopes in influenza
============================================

.. contents::

Summary
----------

This is an analysis of T-cell epitopes in influenza A M1 and NP. Performed by Heather Machkovech in the `Bloom lab`_.

Software used
---------------
* `Python`_ version 2.7.8

* `Biopython`_ version 1.58

* `EMBOSS needle`_ version 6.6.0

* `RAxML`_ version 8.1.2

* `Path-O-Gen`_ version 1.4

* `epitopefinder`_ version 0.11

Main script
--------------
* ``mainepitopecomparison.py`` : main `Python`_ script runs all steps in the analysis through subscripts listed below. 
A more detailed description of the scrips can be found in each subscript.

Subscripts
-------------
* ``parse_sequences.py``  
* ``select_sequences.py`` 
* ``sbatchepitopes.py``
* ``plotepitopesperprotein.py``
* ``CreateSeparateXML.py`` 
* ``sbatchseparatephylogenies.py`` 
* ``MakeMaxCladeCred.py`` 
* ``raxmlhumanswinecombtree.py`` 
* ``CreateXML.py``
* ``sbatchphylogenies.py``
* ``calcaveragemutationpersite.py`` 
* ``calcavemutationantigenicnonantigenic.py``  
* ``plotavemutationepitopetononepitope.py`` 
* ``simplifytree.py`` 
* ``dndsanalysis.py`` 
* ``fanalysis.py`` 
* ``fplots.py``

Steps in analysis
-------------------

Analysis can be run in the terminal by typing: python mainepitopecomparison.py

1. For each protein segment and host, a group of sequences is selected and aligned that will be used to build phylogenies by ``select_sequences.py`` (which in turn runs ``parse_sequences.py``)

2. The location of epitopes in the protein sequences is determined by running epitopefinder. This is accomplished by ``sbatchepitopes.py``.

3. We summarize the findings of epitope finder by plotting the number of epitopes in each protein and the distribution of epitopes in each protein by ``plotepitopesperprotein.py``. We restrict the rest of the analysis to M1 and NP because they contain by far the most epitopes.

4. We build separate protein phylogenies for the human and swine M1 and NP using BEAST. This is accomplished by ``CreateSeparateXML.py`` (which creates the BEAST input file), ``sbatchseparatephylogenies.py`` (which runs BEAST), and ``MakeMaxCladeCred.py`` (which makes a maximum clade credibility tree).

5. We also build combined protein phylogenies for human and swine M1 and NP to determine the time to common ancestor. This is accomplished by ``raxmlhumanswinecombtree.py`` (makes quick tree to check for anamalous sequences), ``CreateXML.py`` (makes input file for BEAST) and ``sbatchphylogenies.py`` (runs BEAST).

6. We next calculate the average number of substitutions that occur at each amino acid site in M1 and NP. This will be used for later analysis. This is accomplished by ``calcaveragemutationpersite.py``. 

7. We compare the rate of epitope substitution to nonepitope substitution for M1 and NP. A summary plot that contains all the values of this epitope to nonepitopesubstitution rate is also constructed. The calculations are completed by  ``calcavemutationantigenicnonantigenic.py`` and the plotting by
``plotavemutationepitopetononepitope.py``.

8. To test if positive selection in M1 and NP epitopes, we perform dN/dS analysis using datamonkey. We plot our analysis looking at the proportion of epitope vs nonepitope sites with dN/dS >1. There is also a cumulative density plot of the dN/dS values for epitope vs nonepitope sites. This analysis is completed by ``simplifytree.py``, which simplifies the maximum clade credibility tree used as input to datamonkey and``dndsanalysis.py``, which creates the plots.

9. We develop another measure to detect positive selection by calculating how many epitopes are contained in an average substitution (called f). We do this for human and swine for the whole tree and just the trunk. We evaluate if human influenza has more epitopes changing than swine in a substitution by computing fhuman/fswine. We also assess if the trunk has an enrichment of epitope-altering substitutions in comparison to the tree by computing ftrunk/ftree. We evaluate significance of these comparisons by making null distributions and calculating p-values. This is done by ``fanalysis.py`` (calculates f values and randomized f values) and ``fplots.py`` (makes summary plots of f values, fhuman/fswine, and ftrunk/ftree and records p-values). 



.. _`Neumann et al 2009`: http://www.nature.com/nature/journal/v459/n7249/full/nature08157.html
.. _`Influenza Virus Resource`: http://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html
.. _`RAxML`: http://sco.h-its.org/exelixis/web/software/raxml/
.. _`Path-O-Gen`: http://tree.bio.ed.ac.uk/software/pathogen/
.. _`Krasnitz et al 2008`: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2519662/
.. _`Python`: https://www.python.org/
.. _`Biopython`: http://biopython.org/wiki/Main_Page
.. _`Path-O-Gen`: http://tree.bio.ed.ac.uk/software/pathogen/
.. _`RAxML`: http://sco.h-its.org/exelixis/web/software/raxml/
.. _`EMBOSS needle`: http://www.ebi.ac.uk/Tools/psa/emboss_needle/
.. _`dos Reis et al 2009`: http://www.ncbi.nlm.nih.gov/pubmed/19787384
.. _`Bloom lab`: http://research.fhcrc.org/bloom/en.html
.. _`epitopefinder`: https://github.com/jbloom/epitopefinder
