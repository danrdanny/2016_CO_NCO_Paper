#!/usr/bin/perl

# This script calculates the number NCO events on the same chromosome arms within
# 4Mb of each other.
# Output is, by column:
# - The number of events observed. Value of 10 means 10 NCO events within 4Mb were observed.
# - The # of trials this number was observed in. A value of 22 means the value was seen in 22 of X trials.
# - The cumulative sum of trials observed.
# - Percent of trials examined. Allows you to determine a confidence interval.

# Danny Miller 

use strict;

my $ncos = 291;
my $stocks = 196;
my $trials = 100000;

my %chrSizes;
$chrSizes{1} = 23540907; 
$chrSizes{2} = 22250000; 
$chrSizes{3} = 25285671; 
$chrSizes{4} = 28109041; 
$chrSizes{5} = 32078544; 

my %countevents;

foreach my $trial (1..$trials) {
	my %output;
	my %count;
	foreach my $gcNum (1..$ncos) {
		my $stock = int(rand($stocks)) + 1;
		my $chrNum = int(rand(5)) + 1;
		my $maxbp = $chrSizes{$chrNum};
		my $pos = int(rand($maxbp)) + 1;

		$output{$stock}{$chrNum} .= "$pos,";
		$count{$chrNum}++;
	}

	my $closeNCOs;
	foreach my $stock (1..$stocks) {
		foreach my $chr (1..5) {
			my $lastNCO;
			foreach my $nco (split /\,/, $output{$stock}{$chr}) {
				if ($lastNCO) {
					my $dist = abs($nco - $lastNCO);
					$closeNCOs++ if $dist < 4000000;
				}
				$lastNCO = $nco;
			}
		}
	}

	$countevents{$closeNCOs}++;
}

print "CloseNCOEvents\tTrialsSeenIn\tCumulativeSum\tPercentageOfEvents\n";
my $sum;
foreach my $key (sort {$a<=>$b} keys %countevents) {
        $sum += $countevents{$key};
        my $percent = sprintf("%0.0f", ($sum / $trials) * 100);
        print "$key\t$countevents{$key}\t\t$sum\t\t$percent\n";
}
