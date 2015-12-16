#!/usr/bin/perl

# This script calculates the number chromsome arms with two or more NCO events. It also
# separately reports chromsome arms with 3 NCO events.
# Output is, by column:
# - The number of events observed. Value of 10 means 10 events matching the criteria were observed.
# - The # of trials this number was observed in. A value of 22 means the value was seen in 22 of X trials.
# - The cumulative sum of trials observed.
# - Percent of trials examined. Allows you to determine a confidence interval.

# Danny Miller

use strict;

my $ncos = 263;
my $stocks = 196;
my $trials = 100000;

my %chrSizes;
$chrSizes{1} = 23540907; 
$chrSizes{2} = 22250000; 
$chrSizes{3} = 25285671; 
$chrSizes{4} = 28109041; 
$chrSizes{5} = 32078544; 

my %countevents;
my %counttriple;
foreach my $trial (1..$trials) {
	my(%output,%count);

	foreach my $gcNum (1..$ncos) {
		my $stock = int(rand($stocks)) + 1;
		my $chrNum = int(rand(5)) + 1;

		$count{$stock}{$chrNum}++;
	}

	my($multiple,$triple) = (0,0);
	foreach my $stock (keys %count) {
		foreach my $chr (keys %{$count{$stock}}) {
			next unless $count{$stock}{$chr} > 1;
			++$multiple;
			$triple++ if $count{$stock}{$chr} == 3;
		}
	}

	$countevents{$multiple}++;
	$counttriple{$triple}++;
}

## Print out summary of chromatids with more than one NCO
print "Count of chromatids with 2 or more NCO events:\n";
print "Chromatids\tTrialsSeenIn\tCumulativeSum\tPercentageOfEvents\n";
my $sum;
foreach my $key (sort keys %countevents) {
	$sum += $countevents{$key};
	my $percent = sprintf("%0.0f", ($sum / $trials) * 100);
	print "$key\t$countevents{$key}\t$sum\t$percent\n";
}

print "\n";

## Print out summary of chromatids with three NCO events
print "Count of chromatids with 3 NCO events:\n";
print "Chromatids\tTrialsSeenIn\tCumulativeSum\tPercentageOfEvents\n";
my $sum;
foreach my $key (sort {$a<=>$b} keys %counttriple) {
	$sum += $counttriple{$key};
	my $percent = sprintf("%0.0f", ($sum / $trials) * 100);
	print "- $key\t$counttriple{$key}\t$sum\t$percent\n";
}
