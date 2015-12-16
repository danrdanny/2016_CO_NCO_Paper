#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options
my %opts;
getopts('a:h', \%opts); # Values in %opts
	
## Usage Statement
if ($opts{'h'} || !$opts{'a'}) {
	print "

	usage: perl co-ncoAnalysis.pl -a <analysis number to run>

	This script returns basic results for the paper.

	Required flags:

	-a  Analysis you want to run (no default). Options are:

	     1 - Per chromosome SNP count
	     2 - Count all events
	     3 - CO events in proximal half vs distal half
	     4 - NCO events in proximal third vs distal two thirds
	     5 - Chromatids with more than one NCO
	     6 - Chromatid with CO and NCO
	     7 - DCO Distribution

	Optional flags: 

	-h  This helpful help.
	\n";
	exit 0;
}

my %chrSizes;
$chrSizes{'chrX'} = 23540907;
$chrSizes{'chr2L'} = 22250000;
$chrSizes{'chr2R'} = 25285671;
$chrSizes{'chr3L'} = 28109041;
$chrSizes{'chr3R'} = 32078544;

my %allStocks;
my $count = 1;
open INF,"./sampleSheet.tsv" or die "Can't open file sampleSheet.tsv: $!";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	next unless $F[4] eq "offspring";

	$allStocks{$F[0]} = $count;
	++$count;
}
close INF;
$count--;
print "[".localtime(time)."] Total stocks identified: $count\n";


## Option 1 - Total and per-chromosome snp counts
if ($opts{'a'} == 1) { 
	my($countVariants,%countSNPs);
	print "[".localtime(time)."] Opening out_uniqueParentalVariants.tsv to return total and per-chromosome SNP counts.\n";
	open INF,"./out_uniqueParentalVariants.tsv" or die "Can't open out_uniqueParentalVariants.tsv: $!";
	while (<INF>) {
        	my(@F) = split /\t/, $_;
        	next unless $F[1] =~ /[0-9]/;

        	my($chr,$id,$ref,$wcns,$wvcf,$cscns,$csvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

		next unless $chrSizes{$chr} > 0;
        	next if $wvcf ne $wcns && $wvcf !~ /\./;
        	next if $csvcf ne $cscns && $csvcf !~ /\./;
        	next if $wcns !~ /[A|G|C|T]/ || $cscns !~ /[A|G|C|T]/;
        	next if $wcns eq $cscns;

        	++$countVariants;
        	$countSNPs{$chr}++;
	}
	print "[".localtime(time)."] $countVariants total variants identified.\n";
	printf "[".localtime(time)."] \n";
	printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
	foreach my $chr (sort keys %chrSizes) {
        	my $snpPerBP = sprintf("%0.0f", $chrSizes{$chr} / $countSNPs{$chr});
        	printf "[".localtime(time)."] %17s %10d - 1 snp / %4d bp\n", $chr, $countSNPs{$chr}, $snpPerBP;
	}

} elsif ($opts{'a'} == 2) {  #Count all events
	my %nonCOChromatids;
	foreach my $stock (keys %allStocks) {
		foreach my $chr (keys %chrSizes) {
			$nonCOChromatids{$stock}{$chr} = 1;
			#print "$stock\t$chr\t$ncos{$stock}{$chr}\n";
		}
	}

	my($ncoCount,$scoCount,$dcoCount,$tcoCount);

	# Count Events
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;

		$scoCount++ if $F[3] eq "sco";
		$dcoCount++ if $F[3] eq "dco";
		$tcoCount++ if $F[3] eq "tco";
		$ncoCount++ if $F[3] eq "nco";
		
		$nonCOChromatids{$F[1]}{$F[0]} = 2 unless $F[3] eq "nco";
	}
	close INF;
	my $totalCount = $scoCount + $dcoCount + $tcoCount;

	my $nonCOCount;
	foreach my $stock (keys %nonCOChromatids) {
		foreach my $chr (keys %chrSizes) {
			$nonCOCount++ if $nonCOChromatids{$stock}{$chr} == 1;
		}
	}

	print "\n";
	print "Non-CO chromatids: $nonCOCount\n";
	print "\n";
	print "SCO: $scoCount\n";
	print "DCO: $dcoCount (" . $dcoCount / 2 ." events)\n";
	print "TCO: $tcoCount (" . $tcoCount / 3 ." events)\n";
	print "Total COs: $totalCount\n";
	print "SCOs + DCO events + TCO events + non-CO chromatids should equal 980\n";
	print "\n";
	print "NCO gene conversions: $ncoCount (counts three discontionus tracts as two event each)\n";

## CO events in proximal half vs distal half
} elsif ($opts{'a'} == 3) {
	my %counts;
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next if $F[3] eq "nco";
		next unless $F[3] eq "sco";

		my $chr = $F[0];
		my $halfChr = $chrSizes{$chr} / 2;

		if ($chr =~ /(chrX|chr2L|chr3L)/) {
			$chr .= "proximal" if $F[2] > $halfChr;
			$chr .= "distal" if $F[2] < $halfChr;
		} elsif ($chr =~ /(chr2R|chr3R)/) {
			$chr .= "proximal" if $F[2] < $halfChr;
			$chr .= "distal" if $F[2] > $halfChr;
		}
		
		$counts{$chr}++;
	}
	close INF;

	foreach my $chr (keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);

		print "$chr proximal: $counts{$proximal} ($proxPercent\%), distal: $counts{$distal} ($distalPercent\%)\n";
	}

## DCO Distribution
} elsif ($opts{'a'} == 4) { 
	my %perChrCounts;
	my $trials = 100000;

	foreach my $trial (1..$trials) {
		my %counts;
		foreach my $event (1..52) {
			my $pos = int(rand(100)) + 1;
			   
			my $chr = 0;
			$chr = 1 if $pos >= 0 && $pos <= 18;
			$chr = 2 if $pos > 19 && $pos <= 35;
			$chr = 3 if $pos > 36 && $pos <= 54;
			$chr = 4 if $pos > 55 && $pos <= 76;
			$chr = 5 if $pos > 77 && $pos <= 100;

			$counts{$chr}++;
		}

		foreach my $chr (keys %counts) {
			$perChrCounts{$chr}{$counts{$chr}}++;
		}
	}

	foreach my $chr (1..5) {
		print "Chromosome: $chr\n";
		my $sum;
		foreach my $key (sort {$a<=>$b} keys %{$perChrCounts{$chr}}) {
			$sum += $perChrCounts{$chr}{$key};
			my $percent = sprintf("%0.1f", ($sum / $trials) * 100);
			print "- $key\t$perChrCounts{$chr}{$key}\t$sum\t$percent\n";
		}
	}

## NCO events in proximal third vs distal two thirds
} elsif ($opts{'a'} == 5) { 
	my %counts;
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next unless $F[3] eq "nco";

		my $chr = $F[0];
		my $thirdOfChr = $chrSizes{$chr} * .33;

		if ($chr =~ /(chrX|chr2L|chr3L)/) {
			my $delim = $thirdOfChr * 2; 
			$chr .= "proximal" if $F[2] > $delim;
			$chr .= "distal" if $F[2] < $delim;
		} elsif ($chr =~ /(chr2R|chr3R)/) {
			$chr .= "proximal" if $F[2] < $thirdOfChr;
			$chr .= "distal" if $F[2] > $thirdOfChr;
		}
		
		$counts{$chr}++;
	}
	close INF;

	print "NCOs in proximal 1/3 vs distal 2/3 of chromosome arm:\n";
	foreach my $chr (keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);

		print "$chr\t$total\tproximal:$counts{$proximal} ($proxPercent\%), distal: $counts{$distal} ($distalPercent\%)\n";
	}

	my %counts;
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next if $F[3] eq "nco";
		next unless $F[3] eq "sco";

		my $chr = $F[0];
		my $thirdOfChr = $chrSizes{$chr} * .33;

		if ($chr =~ /(chrX|chr2L|chr3L)/) {
			my $delim = $thirdOfChr * 2; 
			$chr .= "proximal" if $F[2] > $delim;
			$chr .= "distal" if $F[2] < $delim;
		} elsif ($chr =~ /(chr2R|chr3R)/) {
			$chr .= "proximal" if $F[2] < $thirdOfChr;
			$chr .= "distal" if $F[2] > $thirdOfChr;
		}
		
		$counts{$chr}++;
	}
	close INF;

	print "\n";
	print "SCOs in proximal 1/3 vs distal 2/3 of chromosome arm:\n";
	foreach my $chr (keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);
		print "$chr\t$total\tproximal:$counts{$proximal} ($proxPercent\%), distal: $counts{$distal} ($distalPercent\%)\n";
	}

## Chromatids with more than one NCO
} elsif ($opts{'a'} == 6) { 
	my(%counts,%events);
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next unless $F[3] eq "nco";
		#next unless $F[3] eq "sco";

		next if $F[1] eq "cs12.16" && $F[2] == 12698579;
		next if $F[1] eq "cs12.16" && $F[2] == 12698166;
		next if $F[1] eq "w15.2" && $F[2] == 14707395;
		next if $F[1] eq "w4.2" && $F[2] == 23358042;

		my $chr = $F[0];
	
		$counts{$F[1]}{$F[0]}++;
		$events{$F[1]}{$F[0]} .= "$F[2],";
	}
	close INF;

	my($ncoCounts,$closeNCOs);
	foreach my $stock (keys %counts) {
		foreach my $chr (keys %{$counts{$stock}}) {
			my $count = $counts{$stock}{$chr};
			next unless $count > 1;

			#print "$events{$stock}{$chr}\n";

			my $lastNCO;
			foreach my $nco (split /\,/, $events{$stock}{$chr}) {
				if ($lastNCO) {
					my $dist = abs($nco - $lastNCO);
					$closeNCOs++ if $dist < 4000000;
				}
				$lastNCO = $nco;
			}

			$ncoCounts++;
		}
	}
	print "chromatids with 2 or more NCOs: $ncoCounts.\n";

	print "NCO events <4Mb from each other: $closeNCOs.\n";

## Chromatid with CO and NCO
} elsif ($opts{'a'} == 7) { 
	my(%counts,%events);
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		#next unless $F[3] eq "nco";
		#next unless $F[3] eq "sco";

		next if $F[1] eq "cs12.16" && $F[2] == 12698579;
		next if $F[1] eq "w15.2" && $F[2] == 14707395;
		next if $F[1] eq "w4.2" && $F[2] == 23358042;

		my $chr = $F[0];
	
		$events{$F[1]}{$F[0]} .= "$F[3]|$F[2],";
	}
	close INF;

	my($ncoCounts,$closeNCOs,$countEventsOnSameChromatid,$eventCount,$runAve);
	foreach my $stock (keys %events) {
		foreach my $chr (keys %{$events{$stock}}) {
			my $events = $events{$stock}{$chr};
			my $skip = 1;
			$skip = 0 if $events =~ /nco/ && $events =~ /sco/;
			$skip = 0 if $events =~ /nco/ && $events =~ /dco/;
			$skip = 0 if $events =~ /nco/ && $events =~ /tco/;
			next if $skip == 1;

			++$countEventsOnSameChromatid;

			my %specificTypes;
			foreach my $nco (split /\,/, $events) {
				my($type,$id) = $nco =~ /(\w+)\|(\d+)/;
				next unless $id > 0;
				$specificTypes{$type}{$id} = 1;
			}

			my $closeEvent = 0;
			foreach my $ncoid (keys %{$specificTypes{'nco'}}) {
				foreach my $type ("sco","dco","tco") {
					foreach my $eventid (keys %{$specificTypes{$type}}) {
						my $dist = abs($eventid - $ncoid);

						$closeEvent++ if $dist < 4000000;
						++$eventCount;
						$runAve += $dist;
					}
				}
			}
			$closeNCOs++ if $closeEvent > 0;
		}
	}

	my $aveDistBetween = sprintf("%0.0f", $runAve / $eventCount);

	print "\n";
	print "Chromatids with at least one nco and CO event: $countEventsOnSameChromatid\n";
	print "Chromatids with a NCO event within 4Mb of a CO: $closeNCOs\n";
	print "Average distance between NCO and CO event: $aveDistBetween\n";
}


