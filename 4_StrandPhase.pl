#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use List::Pairwise qw'pair';
use Statistics::Basic qw(:all);

####################################################################################################
##                                                                                                ##
## This script takes as an input partial single cell haplotypes of a single chromosomes and phase ##
## them into two consensus haplotypes representing inherited parental genomes.                    ##
##                                                                                                ##
####################################################################################################

##usage##
if(@ARGV < 1) {
	print "\nNo argument submitted: try 4_StrandPhase.pl -h\|-help\n\n";
	exit;
} elsif ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
}


my @files = ();

foreach (@ARGV) { push @files, $_; }

my $count = 0;
my $prev_count = 0;
my $file1;
my $file2;

my %HoA1 = ();
my %HoA2 = ();

#For a given chromosome find pair of haplotypes with most overlaps and with highest number of mismatched heterozygous positions 
warn "Initializing consensus Haplotypes...\n";
foreach my $a (0..$#files) {
    foreach my $b ($a+1..$#files) {

	open my ($fh1), '<', $files[$a] or die "Can't open file: $!\n";
	open my ($fh2), '<', $files[$b] or die "Can't open file: $!\n";

	my %hap1 = map { (split /\t/)[0,4] } <$fh1>;
	my %hap2 = map { (split /\t/)[0,4] } <$fh2>;

	my @intersection = grep { exists $hap2{$_} } keys %hap1;

	foreach my $pos (@intersection) {
		chomp $hap1{$pos};
		chomp $hap2{$pos};
		
		$count++ if $hap1{$pos} ne $hap2{$pos};

	}
	if ($count > $prev_count) {
		$file1 = $files[$a];
		$file2 = $files[$b];
		$prev_count = $count;
	}
	$count = 0;
	close $fh1;
	close $fh2;
    }
}
warn "DONE...\n";

#haplotypes with most overlaps and with highest number of mismatched heterozygous positions are used as an anchor haplotypes to initialize consensus haplotypes
open my ($fh1), '<', $file1 or die "Can't open file: $!\n";
open my ($fh2), '<', $file2 or die "Can't open file: $!\n";

my %hap1 = ();
my %hap2 = ();
while (<$fh1>) {
	chomp;
	my ($pos, $cons_base, $mapq, $baseq) = (split /\t/)[0,4,5,6];
	my @mapq = (split ",", $mapq);
	my @baseq = (split ",", $baseq);
	my $cov = scalar(@baseq);
	my $bases = ($cons_base x $cov);
	my @bases = (split "", $bases); 		

	push @{$hap1{$pos}{'bases'}}, @bases;
	push @{$hap1{$pos}{'mapq'}}, @mapq;
	push @{$hap1{$pos}{'baseq'}}, @baseq;
	push @{$hap1{$pos}{'files'}}, $file1;
}

while (<$fh2>) {
	chomp;
	my ($pos, $cons_base, $mapq, $baseq) = (split /\t/)[0,4,5,6];
	my @mapq = (split ",", $mapq);
	my @baseq = (split ",", $baseq);
	my $cov = scalar(@baseq);
	my $bases = ($cons_base x $cov);
	my @bases = (split "", $bases); 		

	push @{$hap2{$pos}{'bases'}}, @bases;
	push @{$hap2{$pos}{'mapq'}}, @mapq;
	push @{$hap2{$pos}{'baseq'}}, @baseq;
	push @{$hap2{$pos}{'files'}}, $file2;
}

close $fh1;
close $fh2;

#delete files used to initialize consensus haplotypes from the files list
my @del_indexes = reverse(grep { $files[$_] eq $file1 or $files[$_] eq $file2 } 0..$#files);

foreach my $item (@del_indexes) {
   splice(@files,$item,1);
}

my %dist1 = ();
my %dist2 = ();
my $highest_val = '';

my @sorted_files = sort(@files);
my @process_later = ();


#go through the rest of the haplotypes and assign them to the correct consensus haplotype
while(scalar(@sorted_files) != 0) {
	#continue with next haplotype with most overlapping SNV position either with consensus haplotype 1 or 2
	foreach my $file (@sorted_files) {  
		open my ($fh), '<', $file or die "Can't open file: $!\n";

		my %file = map { (split /\t/)[0,4] } <$fh>;
		close $fh;
	
		my @intersection1 = grep { exists $hap1{$_} } keys %file;
		my @intersection2 = grep { exists $hap2{$_} } keys %file;

		my $dist1 = scalar @intersection1;
		my $dist2 = scalar @intersection2;
	
	$dist1{$file} = $dist1;
	$dist2{$file} = $dist2;
	}

	my $highest_val1 = (sort { $dist1{$b} <=> $dist1{$a} } keys %dist1)[0];
	my $highest_val2 = (sort { $dist2{$b} <=> $dist2{$a} } keys %dist2)[0];

	if ($dist1{$highest_val1} > $dist2{$highest_val2} ) {
		$highest_val = $highest_val1;
	} else {
		$highest_val = $highest_val2;
	}

	#delete haplotype set as next haplotype with most overlaps 
	my ($del_index) = reverse(grep { $sorted_files[$_] =~ /$highest_val/ } 0..$#sorted_files);	
	splice(@sorted_files,$del_index,1);

	#skip haplotypes with less then 5 overlaps
	next if $dist1{$highest_val1} < 5 and $dist2{$highest_val2} < 5;
	%dist1 = ();
	%dist2 = ();

	open my ($fh), '<', $highest_val or die "Can't open file: $!\n";

	my %file = ();
	while (<$fh>) {
		chomp;
		my ($pos, $cons_base, $mapq, $baseq) = (split /\t/)[0,4,5,6];
		my @mapq = (split ",", $mapq);
		my @baseq = (split ",", $baseq);
		my $cov = scalar(@baseq);
		my $bases = ($cons_base x $cov);
		my @bases = (split "", $bases); 		

		push @{$file{$pos}{'bases'}}, @bases;
		push @{$file{$pos}{'mapq'}}, @mapq;
		push @{$file{$pos}{'baseq'}}, @baseq;
		push @{$file{$pos}{'filename'}}, $highest_val;
	}
	close $fh;

	my @intersection1 = grep { exists $hap1{$_} } keys %file;
	my @intersection2 = grep { exists $hap2{$_} } keys %file;

	my $match1 = 0;
	my $match2 = 0;
	my $mismatch1 = 0;
	my $mismatch2 = 0;
	my $count1 = 0;
	my $count2 = 0;

	foreach my $pos (@intersection1) {
		my @bases = @{$file{$pos}{'bases'}};
		my @mapq = @{$file{$pos}{'mapq'}};
		my @baseq = @{$file{$pos}{'baseq'}};
		my $mean_mapq = mean(@mapq);
		my $mean_baseq = mean(@baseq);
		#next if $mean_mapq < 30 or $mean_baseq < 40;
		my $curr_base = '';
		if (@bases == grep { $_ eq $bases[0] } @bases) {
			$curr_base = $bases[0];
		} else { next; }	

		my $cons_base_hap1 = '';
		my @bases_hap1 = @{$hap1{$pos}{'bases'}};
		my @mapq_hap1 = @{$hap1{$pos}{'mapq'}};
		my @baseq_hap1 = @{$hap1{$pos}{'baseq'}};
		my $mean_mapq_hap1 = mean(@mapq_hap1);
		my $mean_baseq_hap1 = mean(@baseq_hap1);
		#next if $mean_mapq_hap1 < 30 or $mean_baseq_hap1 < 40;
		if (@bases_hap1 == grep { $_ eq $bases_hap1[0] } @bases_hap1) {
			$cons_base_hap1 = $bases_hap1[0];
		} else { next; }
		
		$count1++;
		
		if ($curr_base ne $cons_base_hap1) {
			$mismatch1++;
		} else {
			$match1++;
		}

	}

	foreach my $pos (@intersection2) {
		my @bases = @{$file{$pos}{'bases'}};
		my @mapq = @{$file{$pos}{'mapq'}};
		my @baseq = @{$file{$pos}{'baseq'}};
		my $mean_mapq = mean(@mapq);
		my $mean_baseq = mean(@baseq);
		#next if $mean_mapq < 30 or $mean_baseq < 40;
		my $curr_base = '';
		if (@bases == grep { $_ eq $bases[0] } @bases) {
			$curr_base = $bases[0];
		} else { next; }

		my $cons_base_hap2 = '';
		my @bases_hap2 = @{$hap2{$pos}{'bases'}};
		my @mapq_hap2 = @{$hap2{$pos}{'mapq'}};
		my @baseq_hap2 = @{$hap2{$pos}{'baseq'}};
		my $mean_mapq_hap2 = mean(@mapq_hap2);
		my $mean_baseq_hap2 = mean(@baseq_hap2);
		#next if $mean_mapq_hap2 < 30 or $mean_baseq_hap2 < 40;
		if (@bases_hap2 == grep { $_ eq $bases_hap2[0] } @bases_hap2) {
			$cons_base_hap2 = $bases_hap2[0];
		} else { next; }

		$count2++;

		if ($curr_base ne $cons_base_hap2) {
			$mismatch2++;
		} else {
			$match2++;
		}
	}

	next if $count1 == 0 or $count2 == 0;
	
	my $diff1 = ($mismatch1/$count1)*100;
	my $diff2 = ($mismatch2/$count2)*100;
	my $diff = abs($diff1-$diff2);

	if ($diff1 == 0 and $diff2 == 0) {
		push @process_later, $highest_val;
		next;
	}

	my $perc_diff = ((abs($diff1-$diff2))/($diff1+$diff2)/2)*100;
	
	warn "Processing $highest_val...\n";	
	#print "$highest_val\n";
	#print "$match1 $mismatch1 $count1 $diff1\n";
	#print "$match2 $mismatch2 $count2 $diff2\n";
	#print "$perc_diff\n";

	next if $diff1 == $diff2 and push @process_later, $highest_val;
	next if $perc_diff < 25 and push @process_later, $highest_val;
	next if ($mismatch1 + $mismatch2) < 4 and push @process_later, $highest_val;
	

	if ($diff1 < $diff2) { 
		foreach my $pos (keys %file) {
			my @bases = @{$file{$pos}{'bases'}};
			my @mapq = @{$file{$pos}{'mapq'}};
			my @baseq = @{$file{$pos}{'baseq'}};
			my ($filename) = @{$file{$pos}{'filename'}};
				
			push @{$hap1{$pos}{'bases'}}, @bases;
			push @{$hap1{$pos}{'mapq'}}, @mapq;
			push @{$hap1{$pos}{'baseq'}}, @baseq;
			push @{$hap1{$pos}{'files'}}, $filename;
		}

	} else {
		foreach my $pos (keys %file) {
			my @bases = @{$file{$pos}{'bases'}};
			my @mapq = @{$file{$pos}{'mapq'}};
			my @baseq = @{$file{$pos}{'baseq'}};
			my ($filename) = @{$file{$pos}{'filename'}};
				
			push @{$hap2{$pos}{'bases'}}, @bases;
			push @{$hap2{$pos}{'mapq'}}, @mapq;
			push @{$hap2{$pos}{'baseq'}}, @baseq;
			push @{$hap2{$pos}{'files'}}, $filename;
		}
	}

$match1 = 0;
$match2 = 0;
$mismatch1 = 0;
$mismatch2 = 0;
$count1 = 0;
$count2 = 0;

}
	

open my ($out1), '>', 'haplo_1' or die "Can't open file: $!\n";
open my ($out2), '>', 'haplo_2' or die "Can't open file: $!\n";

print $out1 "Pos\tConsensus_base\tCoverage\tBases\tMapq\tBaseq\tCells_IDs\tWeigth\tEntropy\n";
print $out2 "Pos\tConsensus_base\tCoverage\tBases\tMapq\tBaseq\tCells_IDs\tWeigth\tEntropy\n";

open my ($out3), '>', 'Phasing_info_perCell' or die "Can't open file: $!\n";

my %hap1_files = ();
foreach my $pos (sort {$a <=> $b} keys %hap1) {
	my @bases = @{$hap1{$pos}{'bases'}};
	my @mapq = @{$hap1{$pos}{'mapq'}};
	my @baseq = @{$hap1{$pos}{'baseq'}};
	my @files = @{$hap1{$pos}{'files'}};
	
	my $weight = 0;
	my $cons_base = '';
	my %bases = ();	
	my $cov = 0;			
		
	%bases = map {$_, ++$bases{$_}} @bases;
	$cons_base = (sort { $bases{$b} <=> $bases{$a} } keys %bases)[0];
	$cov = $bases{$cons_base};
	if ($cov > 1) { $weight += $cov -1; }

	foreach my $mapq (@mapq) {
		$weight++ if $mapq > 30;
	}
	
	foreach my $baseq (@baseq) {
		$weight++ if $baseq > 40;
	}

	if (scalar(@files) > 1) { $weight += scalar(@files) - 1; }
	
	my $ent_1 = entropy(@bases);
	
	print $out1 "$pos\t$cons_base\t$cov\t@bases\t@mapq\t@baseq\t@files\t$weight\t$ent_1\n";

	foreach my $file (@files) {
		$hap1_files{$file}++;
	}
}

foreach my $filename (sort keys %hap1_files) {
	my ($part1, $part2, $chr, $strand) = (split "_", $filename)[0,1,3,4];
	my $file = join "_", ($part1, $part2);
	print $out3 "$chr\t$file.bam\t";
	print $out3 "w\t" if $strand eq 'A';
	print $out3 "c\t" if $strand eq 'B';
	print $out3 "hap1\n";	 
}

my %hap2_files = ();
foreach my $pos (sort {$a <=> $b} keys %hap2) {
	my @bases = @{$hap2{$pos}{'bases'}};
	my @mapq = @{$hap2{$pos}{'mapq'}};
	my @baseq = @{$hap2{$pos}{'baseq'}};
	my @files = @{$hap2{$pos}{'files'}};
	
	my $weight = 0;
	my $cons_base = '';
	my %bases = ();	
	my $cov = 0;

	%bases = map {$_, ++$bases{$_}} @bases;
	$cons_base = (sort { $bases{$b} <=> $bases{$a} } keys %bases)[0];
	$cov = $bases{$cons_base};
	if ($cov > 1) { $weight += $cov -1; }

	foreach my $mapq (@mapq) {
		$weight++ if $mapq > 30;
	}
	
	foreach my $baseq (@baseq) {
		$weight++ if $baseq > 40;
	}

	if (scalar(@files) > 1) { $weight += scalar(@files) - 1; }
	
	my $ent_2 = entropy(@bases);
	
	print $out2 "$pos\t$cons_base\t$cov\t@bases\t@mapq\t@baseq\t@files\t$weight\t$ent_2\n";

	foreach my $file (@files) {
		$hap2_files{$file}++;
	}
}


foreach my $filename (sort keys %hap2_files) {
	my ($part1, $part2, $chr, $strand) = (split "_", $filename)[0,1,3,4];
	my $file = join "_", ($part1, $part2);
	print $out3 "$chr\t$file.bam\t";
	print $out3 "w\t" if $strand eq 'A';
	print $out3 "c\t" if $strand eq 'B';	
	print $out3 "hap2\n"; 
}

warn "Phasing DONE...\n";

if (@process_later) {
	open my ($out4), '>', 'Discordant_cells' or die "Can't open file: $!\n";
	print $out4 "@process_later\n";
} 

######################################################################################
				## Subroutines ##
######################################################################################

sub entropy {

	my @bases = @_;

	my ($A, $C, $G, $T) = (0,0,0,0);
	my $all_bases = scalar(@bases);

	foreach my $base (@bases) {
		$A++ if $base eq 'A';
		$C++ if $base eq 'C';
		$G++ if $base eq 'G';
		$T++ if $base eq 'T';

	}

	my $A_prob = $A/$all_bases;
	my $C_prob = $C/$all_bases;
	my $G_prob = $G/$all_bases;
	my $T_prob = $T/$all_bases;

	my $ent = -( $A_prob*log2($A_prob) + $C_prob*log2($C_prob) + $G_prob*log2($G_prob) + $T_prob*log2($T_prob) );

	sub log2 {
		my $n = shift;
		return 0 if $n == 0;
		return log($n)/log(2);
	}
	return $ent;
}

sub help { 
print "\n4_StrandPhase.pl <file1 file2 file3...>\n\n";
print "Processes partial single cell haplotypes for a single chromosome\n";
print "stored in the current directory\n";
print "Use mask *singleCellHaps.txt\n\n";
}


