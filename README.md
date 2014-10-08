Phylostratigraphic analysis
=================

####Phylostratigraphy####


Phylostratigraphy was introduced by Domazet-Lo&scaron;o et al. in 2007. The aim of this procedure is to group genes by their phylogenetic origin to uncover footprints of impotant adaptive events in evolution.
It is a statistical approach to reconstruct macroevolutionary trends based on the assumption of punctuated emergence of protein families.


==================

This project builds up a pipeline for phylostratigraphy. It uses BLAST to trace the evolutionary origins based on the curated  non-redundant database of <a href="http://www.ncbi.nlm.nih.gov/">NCBI</a>. The pipeline also implements parsers for BLAST output and an interface to <a href="http://www.monetdb.com/">MonetDB</a> for further analysis.
