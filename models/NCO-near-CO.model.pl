#!/usr/bin/perl

# This script calculates the number NCO events on the same chromosome arms within
# 4Mb of a CO event.
# Output is, by column:
# - The number of events observed. Value of 10 means 10 NCO events within 4Mb of a CO were observed.
# - The # of trials this number was observed in. A value of 22 means the value was seen in 22 of X trials.
# - The cumulative sum of trials observed.
# - Percent of trials examined. Allows you to determine a confidence interval.

# Danny Miller

use strict;

my $ncos = 291;
my $cos = 541; 
my $stocks = 196;
my $trials = 100000;

my %chrSizes;
$chrSizes{1} = 23540907; 
$chrSizes{2} = 23513347; 
$chrSizes{3} = 25285671; 
$chrSizes{4} = 28109041; 
$chrSizes{5} = 32078544; 

my(%eventsPerChrCount,%eventCount,$totalAveDistance);

foreach my $trial (1..$trials) {
	my %cos;
	my %ncos;

	## First populate CO events
	foreach my $coNum (1..$cos) {
		my $stock = int(rand($stocks)) + 1;
		my $chrNum = int(rand(5)) + 1;
		my $maxbp = $chrSizes{$chrNum};
		my $pos = int(rand($maxbp)) + 1;

		$cos{$stock}{$chrNum} .= "$pos,";
	}

	## Populate GC events
	foreach my $gcNum (1..$ncos) {
		my $stock = int(rand($stocks)) + 1;
		my $chrNum = int(rand(5)) + 1;
		my $maxbp = $chrSizes{$chrNum};
		my $pos = int(rand($maxbp)) + 1;

		$ncos{$stock}{$chrNum} .= "$pos,";
	}

	my($curreventsClose4MB,$eventCount,$aveDist) = (0,0,0);

	foreach my $stock (sort keys %cos) {
		foreach my $chr (sort keys %{$cos{$stock}}) {
			my @gcEvents = split /\,/, $ncos{$stock}{$chr};
			my @coEvents = split /\,/, $cos{$stock}{$chr};

			next unless @gcEvents > 0 && @coEvents > 0;

			foreach my $GC (@gcEvents) {
				foreach my $CO (@coEvents) {
					my $gap = abs($GC - $CO);

					$curreventsClose4MB++ if $gap < 4000000; 
					$aveDist += $gap;
					$eventCount++;
				}
			}
		}
	}

	$eventCount{$curreventsClose4MB}++;
	my $totalAve = $aveDist / $eventCount;
	$totalAveDistance += $totalAve;
}

print "NCOEventsCloseToCOs\tTrialsSeenIn\tCumulativeSum\tPercentageOfEvents\n";
my $sum;
foreach my $key (sort {$a<=>$b} keys %eventCount) {
        $sum += $eventCount{$key};
        my $percent = sprintf("%0.1f", ($sum / $trials) * 100);
        print "$key\t$eventCount{$key}\t$sum\t$percent\n";
}

$totalAveDistance = sprintf("%0.0f", $totalAveDistance / $trials);
print "Average distance between COs and NCOs: $totalAveDistance bp\n";
