#!/usr/bin/perl

# This script calculates the number non-crossover, SCO, DCO, TCO, QCO, and greater events expected by
# random chance from 196 individuals.
# The script outputs a summary and also records every trial in an output file.

# Danny Miller

use strict;

my $stocks = 196;
my $totalCOs = 541;
my $arms = 5;
my $trials = 100000;
my $totalArms = $stocks * $arms;

my $output = "Trial\tNon-CO_Count\tSCO_Count\tDCO_Count\tTCO_Count\tQCO_Count\tGreater\n";
my $outputFile = "CO-count-per-arm.output.tsv";
my %output;

foreach my $trial (1..$trials) {
	my %cos;
	foreach (1..$totalArms) {
		$cos{$trial} = 0;
	}
	
	foreach (1..$totalCOs) {
		my $arm = int(rand($totalArms));

		$cos{$arm}++;
	}

	my($nco,$sco,$dco,$tco,$qco,$more);
	foreach (1..$totalArms) {
		$nco++ if $cos{$trial} == 0;
		$sco++ if $cos{$trial} == 1;
		$dco++ if $cos{$trial} == 2;
		$tco++ if $cos{$trial} == 3;
		$qco++ if $cos{$trial} == 4;
		$more++ if $cos{$trial} > 4;
	}

	$output .= "$_\t$nco\t$sco\t$dco\t$tco\t$qco\t$more\n"; 
	$output{'nco'} += $nco;
	$output{'sco'} += $sco;
	$output{'dco'} += $dco;
	$output{'tco'} += $tco;
	$output{'qco'} += $qco;
	$output{'more'} += $more;
}

open OUTF,">./$outputFile";
print OUTF $output;
close OUTF;

print "All output in $outputFile. Total counts below:\n";
print "non-crossover\t$output{'nco'}\n";
print "SCO\t$output{'sco'}\n";
print "DCO\t$output{'dco'}\n";
print "TCO\t$output{'tco'}\n";
print "QCO\t$output{'qco'}\n";
print "greater\t$output{'more'}\n";
