#!/usr/bin/perl

# This script calculates the average distance between DCOs when they are distributed randomly. It also
# counts the number of closely spaced DCOs (<2Mb).

# Danny Miller

use strict;

my $gcs = 291;
my $cos = 541; 
my $stocks = 196;
my $trials = 100000;

my %chrSizes;
$chrSizes{1} = 23540907; 
$chrSizes{2} = 23513347; 
$chrSizes{3} = 25285671; 
$chrSizes{4} = 28109041; 
$chrSizes{5} = 32078544; 

my($allCount,$allGap);
my $output = "Trial\tDCO_Count\tAve_Gap\tCount_Gap_Under_2MB\n";
my $outputFile = "DCO-distribution.output.tsv";

foreach my $trial (1..$trials) {
	my %cos;
	my %output;

	## First populate CO events
	foreach my $coNum (1..$cos) {
		my $stock = int(rand($stocks)) + 1;
		my $chrNum = int(rand(5)) + 1;
		my $maxbp = $chrSizes{$chrNum};
		my $pos = int(rand($maxbp)) + 1;

		my $existingEvents = $cos{$stock}{$chrNum};
		my $skip = 1 if $existingEvents =~ /\d\,\d/;

		my $gap;
		$gap = abs($existingEvents - $pos) unless $skip == 1;

		$cos{$stock}{$chrNum} .= "$pos,";
	}

	my($smallGap,$totalGap,$DCOCount);
	foreach my $stock (sort keys %cos) {
		foreach my $chr (sort keys %{$cos{$stock}}) {
			my $eventCount;
			my @coEvents = split /\,/, $cos{$stock}{$chr};

			next unless @coEvents > 1;
			@coEvents = sort {$a<=>$b} @coEvents;

			my $lastpos;
			foreach my $CO (@coEvents) {
				++$eventCount;
				if ($lastpos) {
					my $gap = abs($CO - $lastpos);

					$smallGap++ if $gap < 2000000;
					$totalGap += $gap;
					++$DCOCount;
				}
				$lastpos = $CO;
			}
		}
	}
	my $aveGap = sprintf("%0.0f", ($totalGap/$DCOCount));
	$output .= "$trial\t$DCOCount\t$aveGap\t$smallGap\n";
	$allCount += $DCOCount;
	$allGap += $aveGap;
}

open OUTF,">./$outputFile";
print OUTF $output;
close OUTF;

print "All data in $outputFile\n";
print "Trials: $trials\n";
print "Ave DCOs/Trial: ".$allCount / $trials."\n";
print "Ave Gap/Trial: ".$allGap / $trials."\n";
