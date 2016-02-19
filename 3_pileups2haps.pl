#!/usr/bin/perl -w
use strict;
use Data::Dumper;

######################################################################################################
##                                                                                                  ##
## This script takes as an input samtools mpilep files and for each position which is different     ##
## from the reference will substitute this base with the observed base at any given position.       ##
## Moreover base and mapping qualities are translated from characters into numbers for further use. ##
## Resulting files are stored into folder per chromosome.                                           ##
##                                                                                                  ##
######################################################################################################

##usage##
if(@ARGV < 1) {
	help();
	exit;
} elsif ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
} else {
	print "WRONG argument submitted: try 3_pileups2haps.pl -h\|-help\n";
	exit;
}

my @files = ();
foreach (@ARGV) { push @files, $_; } #take in pileup files for postprocessing

@files = sort @files;

my %HoA = ();
my %HoA2 = ();
my %HoA3 = ();
my $prefix = '';
my $sanger=33; # Sanger quality=33 # solexa/illumina=64 

foreach my $file ( @files ) {
		
	$prefix = $file;
	$prefix =~ s/\.txt$//;

	my $chr = (split "_", $prefix)[3];
	my $curr_dir = "chr$chr";

	if (not -d $curr_dir) {
		mkdir($curr_dir, 0777) || die "Cannot mkdir newdir: $!";
	}
	
	open IN, '<', $file or die "Can't read from file";

	while ( <IN> ) {
		chomp;
		my ($chr,$pos,$ref,$read,$base_qual,$mapp_qual) = (split /\t/, $_)[0,1,2,4,5,6];

		my @base_qual = split(//, $base_qual);
		my @mapp_qual = split(//, $mapp_qual);

		my @base_qual_dec = ();
		my @mapp_qual_dec = ();

		# Transform quality values from ASCII into Sanger format
      		for( my $i = 0; $i < scalar @base_qual; $i++ ){
       			$base_qual[$i] = ord($base_qual[$i]) - $sanger ;
			push @base_qual_dec, $base_qual[$i];

			$mapp_qual[$i] = ord($mapp_qual[$i]) - $sanger ;
			push @mapp_qual_dec, $mapp_qual[$i];
      		}
		
			
		if (!exists $HoA{$pos}{$ref}) { 
			$HoA{$pos}{$ref} = [$read];
		}else{
			push @{$HoA{$pos}{$ref}}, $read;
		}

		push @{ $HoA2{$pos} }, @base_qual_dec;	
		push @{ $HoA3{$pos} }, @mapp_qual_dec;	
 
	}

my $cmd = '';
$cmd = "rm $file";
system($cmd);

open OUT, '>', "$curr_dir/$prefix".'_hap.txt' or die "Can't create file";

foreach my $posKey1 (sort {$a <=> $b} keys %HoA) {
	foreach my $refKey2 ( keys %{ $HoA{$posKey1} } ) {		

		my $count = @{ $HoA{$posKey1}{$refKey2} };

		my $base;

		foreach my $read ( @{ $HoA{$posKey1}{$refKey2} } ) {

			if ($read =~ /([a|c|g|t])/ig) {
				$base = uc($1);
			}else{
				$base = uc($refKey2);
			}

		}		
		next if ( grep {$_ eq ''} @{ $HoA{$posKey1}{$refKey2} } );
				
		print OUT "$posKey1\t$refKey2\t@{ $HoA{$posKey1}{$refKey2} }\t$count\t$base\t";
		print OUT join ",", @{ $HoA2{$posKey1} };
		print OUT "\t";
		print OUT join ",", @{ $HoA3{$posKey1} };
		print OUT "\n";
 
	}
}
%HoA = ();
%HoA2 = ();
%HoA3 = ();
}

###############################################################################

sub help { 
print "\n3_pileups2haps.pl <file1 file2 file3...>\n\n";
print "Processes files in the current directory\n";
print "Use mask *pileups.txt\n\n";
}

