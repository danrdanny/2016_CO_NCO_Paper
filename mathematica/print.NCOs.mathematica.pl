#!/usr/bin/perl

use strict;

## Print NCO list for mathematica

my $count;
print "Full NCO list (NCOs < 10kb and no discontinous repair events):\n";
print "gclist = {";

open INF,"../CO_NCO_Data/ncoDetail.r6.tsv";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	next unless $F[2] > 0;
	next unless $F[3] eq "nco";

	# skip discontious repair events
        next if $F[1] eq "cs12.16" && $F[2] == 12698579;
        next if $F[1] eq "cs12.16" && $F[2] == 12698166;
        next if $F[1] eq "w15.2" && $F[2] == 14707395;
        next if $F[1] eq "w15.2" && $F[2] == 14707896;
        next if $F[1] eq "w4.2" && $F[2] == 23358042;
        next if $F[1] eq "w4.2" && $F[2] == 23359212;

	my $leftGap  = $F[9] - $F[7];
	my $rightGap = abs($F[10] - $F[13]);
	my $maxGap   = $F[13] - $F[7];

	next if $maxGap > 9999;

	print "{$maxGap, $leftGap, $rightGap}, " if $leftGap < $rightGap;
	print "{$maxGap, $rightGap, $leftGap}, " if $leftGap > $rightGap;
	++$count;
}

print "}\n";
print "\n";
print "Count of events in the full NCO list: $count\n";
print "\n";
print "\n";

my $count;
print "Reduced NCO list (NCOs with left and right gaps <= 1000 bp and no discontinous repair events):\n";
print "gclist = {";

open INF,"../CO_NCO_Data/ncoDetail.r6.tsv";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	next unless $F[2] > 0;
	next unless $F[3] eq "nco";

	# skip discontious repair events
        next if $F[1] eq "cs12.16" && $F[2] == 12698579;
        next if $F[1] eq "cs12.16" && $F[2] == 12698166;
        next if $F[1] eq "w15.2" && $F[2] == 14707395;
        next if $F[1] eq "w15.2" && $F[2] == 14707896;
        next if $F[1] eq "w4.2" && $F[2] == 23358042;
        next if $F[1] eq "w4.2" && $F[2] == 23359212;

	my $leftGap  = $F[9] - $F[7];
	my $rightGap = abs($F[10] - $F[13]);
	my $maxGap   = $F[13] - $F[7];

	next if $leftGap  > 1000;
	next if $rightGap > 1000;

	print "{$maxGap, $leftGap, $rightGap}, " if $leftGap < $rightGap;
	print "{$maxGap, $rightGap, $leftGap}, " if $leftGap > $rightGap;
	++$count;
}

print "}\n";
print "\n";
print "Reduced NCO count: $count\n";
print "\n";
