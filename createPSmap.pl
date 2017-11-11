#!/usr/bin/perl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W)   Alexander Gabel
# Copyright (C) Alexander Gabel
#
# createPSmap.pl wrapper script to start BLAST and create PS map with java xml parser

use Cwd;
use File::Basename;
use File::Spec;
use Getopt::Long;
use POSIX qw/ceil/;

my $refOrg;
my $db;
my $psFile;
my $seqOffset = 50;
my $evalue = 1e-5;
my $threads = 1;
my $blastPlus;

GetOptions('help' => \$help,
'organism=s' => \$refOrg,
'database=s' => \$db,
'seqOffset=i' => \$seqOffset,
'evalue=s' => \$evalue,
'threads=i' => \$threads,
'prefix=s' => \$outputPrefix,
'blastPlus' => \$blastPlus);

if(defined $help){
  usage();
}

if( ((!defined $refOrg) or (!defined $db)) or (!defined $outputPrefix)){
  print $refOrg."\n";
  print $db."\n";
  print $outputPrefix."\n";
	usage();
}

my $scriptLocation = File::Basename::dirname(Cwd::abs_path($0));

my $sep = File::Spec->catfile('', '');

# Loading genes from fasta file and store in hash
my %seqs = readFASTA($refOrg);
my @seqArr = sort keys %seqs;

# Begin to BLAST and create for each bulk of sequences a phylostatigraphic map
my @path = split("/",$refOrg);
my $prefix = pop(@path);

$prefix =~ s/\.(.+)//ig;

my @psFiles = ();
my @tableFiles = ();

my $i = 0;

my $xmlLoc = "";
my $xmlDir;

while($i < scalar(@seqArr)){
  
  my $outFile = $outputPrefix."_".$prefix."_".($i+1)."_".($i+$seqOffset).".fa";
  
  if($i + $seqOffset > scalar @seqArr){
  		$outFile = $outputPrefix."_".$prefix."_".($i+1)."_".(scalar @seqArr).".fa";
  }
  
  my $xmlFile = $outFile;
  $xmlFile =~ s/\.fa/.xml/;
  
  if(-e $xmlFile.".tbz"){
    
    print "XML file (".$xmlFile.") already exists. Skipping BLAST procedure...\n";
    
    $command = "tar tf ".$xmlFile.".tbz";
    $xmlLoc = `$command`;
    chomp $xmlLoc;

    $command = "tar xfvj ".$xmlFile.".tbz";
    system($command);

    $xmlFile = $xmlLoc;

    if(! defined $xmlDir){
	($vol, $directories,$f) = File::Spec->splitpath( $xmlLoc );
	@dirs = File::Spec->splitdir( $directories );
	$xmlDir = $dirs[0];
    }	    
    
    $i += $seqOffset;
    
  }else{
    
    my @outIds = ();
    for(my $j = 0; $j < $seqOffset; $j++){
      
      $outIds[$j] = $seqArr[$i];
      
      if($i >= (scalar @seqArr - 1)){
      		last;
      }
      
      $i++;
    }
    
    writeFASTA(\%seqs,$outFile,@outIds);
    
    my $command;
    my $blastCommand;
    
    if(defined $blastPlus){
      $blastCommand = "blastp -query ".$outFile." -db ".$db." -num_threads ".$threads." -out ".$xmlFile." -task blastp -evalue ".$evalue." -max_target_seqs 65200 -matrix BLOSUM62 -outfmt 5";
    }else{
      $blastCommand = "blastall -p blastp -i ".$outFile." -d ".$db." -a ".$threads." -m 7 -o ".$xmlFile." -e ".$evalue." -b 65200 -M BLOSUM62";
    }
    
    `$blastCommand`;
    
    $command = "tar cfvj ".$xmlFile.".tbz ".$xmlFile;
    `$command`;
    
    unlink($outFile);
  }
  
  $command = "java -jar ".$scriptLocation.$sep."ParseXMLtoPS.jar -i ".$xmlFile." -e ".$evalue;
  my $jarOut = `$command`;
  
  my @files = split("\n",$jarOut);
  @files = split(" and ",$files[(scalar (@files)) -1]);
  
  push(@tableFiles,$files[0]);
  push(@psFiles,$files[1]);
  
  # deleting uncompressed xml files
  print "Removing ".$xmlFile." after compressing to ".$xmlFile.".tbz\n";
  unlink($xmlFile);
}

# Summarize all subset PS maps together

my $finalMapFile = $outputPrefix."_final_ps_map.csv";

print "Creating final phylostratigraphic map -> ".$finalMapFile."\n";

my $command = "cat ".join(" ",@psFiles)." | sort -u | grep -v \"PS;query_id\" | sed \'1 s/.*/PS;GeneID\\n&/\' > ".$finalMapFile;
`$command`;

# deleting subset PS maps
unlink(@psFiles);


my $tableArchive = $outputPrefix."_BLAST_PS_tables.tbz";

print "Creating an archive file ".$tableArchive." containing all BLAST results represented in tables.\n";

$command = "tar cfvj ".$tableArchive." ".join(" ",@tableFiles);
`$command`;

# deleting uncompressed table files
unlink(@tableFiles);

if(-d $xmlDir){
  my $curWd = `pwd`;
  chomp $curWd;
  print "\n---> !!! Please remove the temporary directory \'".$curWd.$sep.$xmlDir."\' which was created during untaring previous BLAST results!!! <---\n\n";
}

# print usage of Perl script
sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: perl createPSmap.pl [--organism organism_proteom.fasta] [--prefix prefix_output_files] [--seqOffset 50] [--evalue 1e-5] [--threads 1] [--blastPlus] [--help]\n\n";
  print "-o --organism \t\t FASTA file with amino acis sequences containing the proteom\n";
  print "-d --database \t\t BLAST database\n";
  print "-p --prefix \t\t prefix for output files \n";
  print "-s --seqOffset \t\t size of a package the FASTA file should split in (default: 50)\n";
  print "-e --evalue \t\t E-value for BLAST comparisons and assignment of protein coding genes to their phylostratum (default: 1e-5)\n";
  print "-b --blastPlus \t\t BLAST+ is used for similarity searches otherwise \'blastall -p blastp\' is used\n";
  print "-h --help \t\t print this message\n";
  exit;
}

sub readFASTA{
  
  my $file = shift;
  
  my %seqHash = ();
  
  open(IN, $file) or die $!;
  my $id;
  while(my $line = <IN>){
    
    if($line =~ /^>(.*)\n/ ){
      $id = $1;
      $seqHash{$id} = "";
      
    }else{
      chomp($line);
      $seqHash{$id} .= $line;
    }
  }
  close(IN);
  
  return %seqHash;
}

sub writeFASTA{
  
  my $seqHash = shift;
  my $outFile = shift;
  my @idArr = @_;
  
  # if no array of ids is given write the whole seqHash into fasta file
  
  if((scalar @idArr) == 0){
    @idArr = sort keys %$seqHash;
  }
  
  open(OUT,">",$outFile) or die $!;
  for my $id (@idArr){
    
    my $seq = $seqHash->{$id};
    print OUT ">".$id."\n";
    while($seq =~/.{1,80}/){
      print OUT "$&\n";
      $seq = $';
    }
  }
  
}
