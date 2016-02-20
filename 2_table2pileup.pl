#!/usr/bin/perl -w
use strict;
use Data::Dumper;

####################################################################################################################################
##                                                                                                                                ##
## This script will split localized WC regions in single cell Strand-seq libraries by direction they map to the reference genome ##
##                                                                                                                                ##
####################################################################################################################################

##usage##
if(@ARGV < 1) {
	print "WRONG argument submitted: try 2_table2pileup.pl -h\|-help\n";
	exit;
} elsif ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
}

my $genome = shift; #'/home/daewoooo/WORK/GENOMES/NCBI36/NCBI36.fa';
my $method = shift; #single or paired: for single-end or paired-end reads
my $infile = shift; #table containing coordinates of CW regions in all single cells

open IN, '<', $infile or die "Can't read from file";

my $prev_chr = '';
my $count = 1;

while(<IN>) {
	$_ =~ /(\d+):(\d+):(\d+):(.+).bam/; #read in WC region table => W and C reads can be distinguished based on directionality
	my ($chr, $start, $end, $prefix) = ($1, $2, $3, $4);
	$infile = $prefix.'.bam';
	print "$infile\n";
	my $num = (split /_/, $infile)[1];

	#filter out reads with mapping qulaity less than 10
	open F1, "| samtools view -bS - | samtools mpileup -s -q 10 -l snp_list_$chr.txt -f $genome - | cat >> ".$prefix."_mpileup_".$chr."_A.txt";
	open F2, "| samtools view -bS - | samtools mpileup -s -q 10 -l snp_list_$chr.txt -f $genome - | cat >> ".$prefix."_mpileup_".$chr."_B.txt";

	$chr = 'X' if $chr == 23;
	open F, 'samtools view -h '.$infile.' '.$chr.':'.$start.'-'.$end.' |'; #add chr before chromosome ID in case of different bam header

	while ( <F> ) {
		my $flag = '';
   		if ( m/^\@/ ) {
     			print F1;
     			print F2;
   		} else {
     			$flag = (split /\t/, $_)[1];
			next if $flag == 1024; #filter out duplicate reads

			#split reads based on the direction they map to the reference genome
			if ($method eq 'paired') {
     				print F1 $_ if $flag == 99 or $flag == 147;
     				print F2 $_ if $flag == 83 or $flag == 163; 
			} elsif ($method eq 'single') {
				print F1 $_ if $flag == 16;
				print F2 $_ if $flag == 0;
			} else {
				print "WARRNING: wrong method submitted!\n";
			} 
   		}	
	}		
}

##################################################################################

sub help { 
print "\n2_table2pileup.pl <genome> <method> <WC_regions_table>\n\n";
print "<genome> :Reference genome assembly in fasta\n";
print "<method> :single - single-end reads\n";
print "\t :paired - paired-end read\n";
print "<WC_regions_table> :File containing positions of WC regions. (Format chromosome:start:end:bamfileID)\n\n";
}

