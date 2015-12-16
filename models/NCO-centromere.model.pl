#!/usr/bin/perl

# This script calculates the number NCO events within the proximal 1/3 of each chromosome arm.
#
# Output is per chromosome. By column:
# - The number of events observed. Value of 10 means 10 NCO events within the proximal 1/3 were observed.
# - The # of trials this number was observed in. A value of 22 means the value was seen in 22 of $trials trials.
# - The cumulative sum of trials observed.
# - Percent of trials examined. Allows you to determine a confidence interval.

# Danny Miller

use strict;

my $gcs = 291;
my $stocks = 196;
my $trials = 1000;

my %chrSizes;
$chrSizes{1} = 23540907; 
$chrSizes{2} = 22250000; 
$chrSizes{3} = 25285671; 
$chrSizes{4} = 28109041; 
$chrSizes{5} = 32078544; 

my %chrKey;
$chrKey{1} = "chrX";
$chrKey{2} = "chr2L";
$chrKey{3} = "chr2R";
$chrKey{4} = "chr3L";
$chrKey{5} = "chr3R";

my %countevents;
foreach my $trial (1..$trials) {
	my(%total,%proximal);
	foreach my $gcNum (1..$gcs) {
		my $chr = int(rand(5)) + 1;
		my $maxbp = $chrSizes{$chr};
		my $pos = int(rand($maxbp)) + 1;

		my $thirdOfChr = $chrSizes{$chr} * .33;

		my $loc;
		if ($chr =~ /(1|2|4)/) {
			my $delim = $thirdOfChr * 2;
			$loc = "proximal" if $pos > $delim;
			$loc = "distal" if $pos < $delim;
		} elsif ($chr =~ /(3|5)/) {
			$loc = "proximal" if $pos < $thirdOfChr;
			$loc = "distal" if $pos > $thirdOfChr;
		}

		$total{$chr}++;
		$proximal{$chr}++ if $loc =~ /proximal/;
	}

	foreach my $chr (sort keys %chrSizes) {
		my $proxPercent = sprintf("%0.0f",($proximal{$chr} / $total{$chr}) * 100);
		$countevents{$chr}{$proxPercent}++;
	}
}

print "Per-chromosome output:\n";
foreach my $chr (1..5) {
	my $sum;
	print "Number of NCOs within the proximal 1/3 of $chrKey{$chr}:\n";
	foreach my $key (sort {$a<=>$b} keys %{$countevents{$chr}}) {
	        $sum += $countevents{$chr}{$key};
	        my $percent = sprintf("%0.1f", ($sum / $trials) * 100);
	        print "$chrKey{$chr}\t$key\t$countevents{$chr}{$key}\t$sum\t$percent\n";
	}
	print "\n";
}
