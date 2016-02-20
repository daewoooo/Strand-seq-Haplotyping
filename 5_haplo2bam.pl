#!/usr/bin/perl -w
use strict;
use Data::Dumper;

####################################################################################################################################
##                                                                                                                                ##
## This script will split localized WC regions in single cell Strand-seq libraries per direction they map to the reference genome ##
##                                                                                                                                ##
####################################################################################################################################

##usage##
if(@ARGV < 1) {
	rint "WRONG argument submitted: try 5_haplo2bam.pl -h\|-help\n";
	exit;
} elsif ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
}

my $region_tab = shift or die 'No infile';
my $haplo_files = shift or die 'No infile';
my $chrom = shift or die 'ERROR';

open my ($fh1), '<', $region_tab or die "Can't open file: $!\n"; #read in coordinates of CW regions
open my ($fh2), '<', $haplo_files or die "Can't open file: $!\n"; #read in phasing info for single cell libraries 

my %region_table = ();
while (<$fh1>) {
	chomp;
	my ($chr,$start,$end,$filename) = (split ":", $_);
	next if $chr != $chrom;
	push @{ $region_table{$filename} }, [$start,$end];
}

my $counter = 1;
while (<$fh2>) {
	chomp;
	my ($chr, $filename, $strand, $hap) = (split "\t", $_);

	$chrom = 'X' if $chrom == 23;

	if ( exists($region_table{$filename}) ) {
		my $prefix = $filename;
		$prefix =~ s/\.bam$//;

		foreach my $region ( @{$region_table{$filename}} ) {
			my $start = $region->[0];
			my $end = $region->[1];

			open F, 'samtools view -h '.$filename.' '.$chrom.':'.$start.'-'.$end.' |';
			open F1, '| samtools view -bS - > '.$counter.'_'.$prefix.'_'.$hap.'.bam';
			$counter++;

			while (<F>) {
				my $flag = '';
   				if ( m/^\@/ ) {
					print F1;
				} else {
					$flag = (split /\t/, $_)[1];
					#next if $flag == 1024;
					if ($strand eq 'w') {
     						print F1 $_ if $flag == 99 or $flag == 147 or $flag == 16;
					} elsif ($strand eq 'c') {
						print F1 $_ if $flag == 83 or $flag == 163 or $flag == 0;
					} 
				}
			}
			close F;
			close F1;		
		}
	}
}

## merge created regional bam files into a single bam file for a given chromosome
my $cmd;

$cmd = join( ' ', 'samtools', 'merge',
                  'chr'.$chrom.'_haplo1.bam',
                  '*hap1.bam',
               );
system $cmd;

$cmd = join( ' ', 'rm', '*hap1.bam');
system $cmd;

$cmd = join( ' ', 'samtools', 'merge',
                  'chr'.$chrom.'_haplo2.bam',
                  '*hap2.bam',
               );
system $cmd;

$cmd = join( ' ', 'rm', '*hap2.bam');
system $cmd;


##################################################################################

sub help { 
print "\n5_haplo2bam.pl <WC_regions_table> <Phasing_info_perCell> <chromosome>\n\n";
print "<WC_regions_table> :File containing positions of WC regions. (Format chromosome:start:end:bamfileID)\n";
print "<Phasing_info_perCell> :Phasing info of each single cell deduced by 4_StrandPhase.pl\n"; 
print "\t\t\t(Format chromosome BamFileName directionalty(c\|w) phase(hap1\|hap2)\n";
print "<chromosome> : Chromosome number (1,2,3...) \n\n";
}
