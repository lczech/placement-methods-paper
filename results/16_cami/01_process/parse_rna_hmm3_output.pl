#!/usr/bin/env perl
use strict;
use warnings;
$|++;

unless ($ARGV[0]) {die "Usage: parse_rna_hmm_output.pl <rna_hmm3_output> \n"}

my $output = shift @ARGV;
open (IN, "<$output");
open (OUT, ">$output.parsed");
open (OUT2, ">$output.stats");
my $read_id;
my %reads;
my %hits;
my $new_evalue;
my $old_evalue;
my $row;
while (<IN>){
    unless (/^##/){
	$row = $_;
	chomp $row;
	my @line = split "\t",$_;
	$read_id = $line[0];
	$new_evalue = $line[5];
	$reads{$read_id}++;

###################################################
#CHECK if there is more than one hit for this read#
###################################################
	if ($reads{$read_id} > 1){
	    my @line2 = split "\t", $hits{$read_id};
	    $old_evalue = $line2[5];
#	    print "old hit: ".$hits{$read_id}."\n";
#	    print "old evalue: " . $old_evalue."\n";
#	    print "new evalue: " . $new_evalue."\n";
######################################################
##REPLACE old entry w/ new one if new evalue is lower#
######################################################
	    if ($old_evalue > $new_evalue){
#		print "OLD: ".$hits{$read_id}."\n";
		$hits{$read_id} = $row;
#		print "NEW: ".$hits{$read_id}."\n";
	    }
	}
	else {
		$hits{$read_id} = $row;
	}
    }
}
close IN;

#Count number hits
my $num_hits = keys %reads;
$num_hits = $num_hits -1; #substract 1 b/c of header line

my $bac16S = my $bac23S = my $arc16S = my $arc23S = my $euk18S = my $euk28S = 0;
foreach my $i (keys %hits){
    my @line3 = split "\t", $hits{$i};
#    print $line3[8],"\n";
    $bac16S++ if $line3[8] =~ /bac 16S_rRNA/;
    $bac23S++ if $line3[8] =~ /bac 23S_rRNA/;
    $arc16S++ if $line3[8] =~ /arc 16S_rRNA/;
    $arc23S++ if $line3[8] =~ /arc 23S_rRNA/;
    $euk18S++ if $line3[8] =~ /euk 18S_rRNA/;
    $euk28S++ if $line3[8] =~ /euk 28S_rRNA/;
	print OUT $hits{$i}."\n";
}
print OUT2 "sample" . "\t" . $output . "\n";
print OUT2 "# of hits" . "\t" . $num_hits . "\n";
print OUT2 "bac 16S" . "\t" . $bac16S . "\n";
print OUT2 "bac 23S" . "\t" . $bac23S . "\n";
print OUT2 "arc 16S" . "\t" . $arc16S . "\n";
print OUT2 "arc 23S" . "\t" . $arc23S . "\n";
print OUT2 "euk 18S" . "\t" . $euk18S . "\n";
print OUT2 "euk 28S" . "\t" . $euk28S . "\n";
