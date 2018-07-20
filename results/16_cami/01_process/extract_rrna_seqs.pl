#!/usr/bin/env perl
use warnings;
use strict;
$|++;

#USAGE: extract_rrna_seqs.pl <*.rRNA.parsed> <number of splits> <min length>
my $parsed_file = shift @ARGV;
my $split = shift @ARGV; #into how many files are fasta files split?
my $min = shift @ARGV;
$parsed_file =~ m/^(.+).rRNA.parsed/;
my $file = $1;

system "rm $file*_rRNA";

#get read ids from rna_hmm3 parsed output, fetch sequences and reverse complement if necessary
my @rrna = ("bac 16S_rRNA","arc 16S_rRNA","euk 18S_rRNA");

#uncomment next line if lsu should be included
#my @rrna = ("bac 16S_rRNA","bac 23S_rRNA","arc 16S_rRNA","arc 23S_rRNA","euk 18S_rRNA","euk 28S_rRNA");

foreach my $j (@rrna){
  print $j . "\n";
  $j =~ m/(.+) (.+)/;
  my $type = $2;
  open IN,"<$parsed_file";
  open OUT,">>$file.$type";
  while (<IN>){
		my @line = split "\t",$_;
		if (/$j/){
	    my $start = $line[3];
	    my $end = $line[4];
	    my $offset = $start - 1;
#	    my $offset = $start + 24;
	    my $length = $end - $offset;
#	    my $length = $end - $offset - 25;
	    unless ($length < $min){
				for (my $i=0; $i<$split; $i++){
	    		my $command;
	    		if ($split == 1){
			# print $file;<stdin>;
			# print $line[0];<stdin>;
						$command = "cdbyank ./$file.fasta.cidx -d ./$file.fasta -a \"$line[0]\"";
	    		}else{
						$command = "cdbyank ./$file.$i.fasta.cidx -d ./$file.$i.fasta -a \"$line[0]\"";
	    		}
	    		# print $command;<stdin>;
	    		my @sequence = split "\n",`$command`;
	    		unless (scalar @sequence == 0){
						my $name=''; my $seq='';
						foreach my $i (@sequence){
		    			if ($i =~ /^>/){
								$name = $i;
		    			}else{
								chomp $i;
								$seq .= $i;
		    			}
						}
						if ($line[6] eq "-"){
					    chomp($seq);
					    my $string = substr $seq, $offset, $length;
					    my $rev = reverse $string;
					    $rev =~ tr/GATCgatc/CTAGctag/;
					    print OUT $name . "_rc" . "\n" . $rev . "\n";
						}else{
					    my $string = substr $seq, $offset, $length;
					    print OUT $name . "\n" . $string . "\n";
						}
	    		}
				}
    	}
		}
  }
  close IN;
  close OUT;
}
