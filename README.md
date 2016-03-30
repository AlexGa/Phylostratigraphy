Phylostratigraphic analysis
=================

This project builds up a pipeline for phylostratigraphy. It uses BLAST to trace the evolutionary origins based on the curated  non-redundant database of <a href="http://www.ncbi.nlm.nih.gov/">NCBI</a>.

__Phylostratigraphy__ was introduced by <a href="http://www.sciencedirect.com/science/article/pii/S0168952507002995">Domazet-Lo&scaron;o et al. in 2007</a> to trace the evolutionary origin of protein coding genes. Thus, it groups genes by their phylogenetic origin to uncover footprints of important adaptive events in evolution.
It is a statistical approach to reconstruct macroevolutionary trends based on the assumption of punctuated emergence of protein families.

##Performing Phylostratigraphy

It can be performed by using the Perl script `createPSmap.pl`. The resulting phylostratigraphic map stores the phylostratum in the first column and the corresponding gene id in the second column.

For creating the phylostratigraphic map the following steps have to be done:

1) Make sure that BLAST is installed on your machine. You can choose between BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/) and BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

2) Download the sequence database <a href="http://msbi.ipb-halle.de/download/phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz">phyloBlastDB_Drost_Gabel_Grosse_Quint.fa</a> used for BLAST searches and unpack it 
`tar xfvj phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz`

3) Format the BLAST sequence database.<br />
If you are using BLAST<br />
`formatdb -p T -i phyloBlastDB_Drost_Gabel_Grosse_Quint.fa`<br />
or, if you are using BLAST+<br />
`makeblastdb -dbtype prot -in phyloBlastDB_Drost_Gabel_Grosse_Quint.fa`<br />

4) Make sure that the header of your FASTA-file of the organism you are interested in (e.g. Athaliana_167_protein.fa) fullfills the following specification:<br />
<code>>GeneID | [organism_name] | [taxonomy]</code><br />
Notice, the taxonomy begins after the node "Cellular organisms" e.g.
```{terminal}
>NP_146894.1 | [Aeropyrum pernix] | [Archaea; Crenarchaeota; Thermoprotei; Desulfurococcales; Desulfurococcaceae; Aeropyrum]
or
>YP_001514406.1 | [Acaryochloris marina MBIC11017] | [Bacteria; Cyanobacteria; Oscillatoriophycideae; Chroococcales; Acaryochloris; Acaryochloris marina]
or
>ATCG00500.1|PACid:19637947 | [Arabidopsis thaliana] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis]
```

5) Use the following command to start the Perl script<br />
```terminal
perl createPSmap.pl [--organism organism_proteom.fasta] [--seqOffset 50] [--evalue 1e-5] [--threads 1]  [--blastPlus] [--help]

Arguments:
-o --organism     FASTA file with amino acis sequences containing the proteom
-d --database     BLAST database
-p --prefix       prefix for output files 
-s --seqOffset    size of a package the FASTA file should split in (default: 50)
-e --evalue       E-value for BLAST comparisons and assignment of protein coding genes to their phylostratum (default: 1e-5)
-b --blastPlus    BLAST+ is used for similarity searches otherwise 'blastall -p blastp' is used
-h --help         print this message
```

E.g. Starting pipeline for *A. thaliana* using BLAST (`blastall -p blastp`) for similarity searches<br />
```terminal
perl createPSmap.pl --organism Athaliana_167_protein_with_new_Header.fa --seqOffset 50  --evalue 1e-5 --threads 60
```
E.g. Starting pipeline for *A. thaliana* using BLAST+ for similarity searches<br />
```terminal
perl createPSmap.pl --organism Athaliana_167_protein_with_new_Header.fa --seqOffset 50  --evalue 1e-5 --threads 60 --blastPlus
```