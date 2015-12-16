#!/usr/bin/perl

# This script generates 
#
# Danny Miller

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

	     1 - Genome-wide and per-chromosome SNP count
	     2 - Count of SCOs, DCOs, TCOs, NCOs, and non-crossover chromatds
	     3 - CO events in proximal half vs distal half of each chromatid
	     4 - NCO events in proximal third vs distal two thirds of each chromatid
	     5 - Chromatids with more than one NCO
	     6 - Chromatid with CO and NCO

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
open INF,"./co-nco-sampleSheet.tsv" or die "Can't open file co-nco-sampleSheet.tsv: $!";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	next unless $F[4] eq "offspring";

	$allStocks{$F[0]} = $count;
	++$count;
}
close INF;
$count--;
print "[".localtime(time)."] Option passed: $opts{'a'}\n";
print "[".localtime(time)."] Total stocks identified: $count\n\n";

  #--------------------------------------------------#
  #  Option 1 - Total and per-chromosome snp counts  #
  #--------------------------------------------------#

if ($opts{'a'} == 1) { 
	my($countVariants,%countSNPs);
	print "[".localtime(time)."] Genome-wide and per-chromosome SNP counts are from co-nco-uniqueParentalVariants.tsv.\n";
	open INF,"./co-nco-uniqueParentalVariants.tsv" or die "Can't open co-nco-uniqueParentalVariants.tsv: $!";
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
	print "[".localtime(time)."] Total SNPs identified: $countVariants\n";
	print "[".localtime(time)."] \n";
	printf "[".localtime(time)."] Per chromosome SNP counts:\n";
	foreach my $chr (sort keys %chrSizes) {
        	my $snpPerBP = sprintf("%0.0f", $chrSizes{$chr} / $countSNPs{$chr});
        	printf "[".localtime(time)."] %6s %8d\t1 snp / %-4d bp\n", $chr, $countSNPs{$chr}, $snpPerBP;
	}

  #--------------------------------------------#
  # Option 2 - Count NCO, SCO, DCO, etc events #
  #--------------------------------------------#

} elsif ($opts{'a'} == 2) {
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
	print "Count of CO events.\n";
	print "\n";
	print "Non-CO chromatids: $nonCOCount\n";
	print "\n";
	print "SCO: $scoCount\n";
	print "DCO: $dcoCount (" . $dcoCount / 2 ." events)\n";
	print "TCO: $tcoCount (" . $tcoCount / 3 ." events)\n";
	print "Total COs: $totalCount\n";
	print "SCOs + DCO events + TCO events + non-CO chromatids should equal 980\n";
	print "\n";
	print "NCO gene conversions: $ncoCount (counts three discontionus tracts as two events each)\n";
	print "\n";

  #-----------------------------------------------------------#
  # Option 3: Count CO events in proximal half vs distal half #
  #-----------------------------------------------------------#

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

	print "Count of CO events in the proximal 1/2 vs distal 1/2 of each chromosome arm.\n";
	print "\n";

	print "Chr\tProximal_Events\t\tDistal_Events\n";
	foreach my $chr (sort keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);

		print "$chr\t$counts{$proximal} ($proxPercent\%)\t\t$counts{$distal} ($distalPercent\%)\n";
	}

  #------------------------------------------------------------------------------------------#
  # Option 4: Count NCO events in proximal third vs distal two thirds of each chromosome arm #
  #------------------------------------------------------------------------------------------#

} elsif ($opts{'a'} == 4) { 
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
	print "\n";
	print "Chr\tEvents\tProximal_Events\t\tDistal_Events\n";
	foreach my $chr (keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);

		print "$chr\t$total\t$counts{$proximal} ($proxPercent\%)\t\t$counts{$distal} ($distalPercent\%)\n";
	}

	my %counts;
	$counts{"chr2Rproximal"} = 0; #there are none, so we set it to zero to make the output pretty
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
	print "\n";
	print "Chr\tEvents\tProximal_Events\t\tDistal_Events\n";
	foreach my $chr (keys %chrSizes) {
		my $proximal = $chr."proximal";
		my $distal = $chr."distal";
		my $total = $counts{$proximal} + $counts{$distal};

		my $proxPercent = sprintf("%0.1f",($counts{$proximal} / $total) * 100);
		my $distalPercent = sprintf("%0.1f", $counts{$distal} / $total * 100);
		print "$chr\t$total\t$counts{$proximal} ($proxPercent\%)\t\t$counts{$distal} ($distalPercent\%)\n";
	}

  #---------------------------------------------------#
  # Option 5: Count chromatids with more than one NCO #
  #---------------------------------------------------#

} elsif ($opts{'a'} == 5) { 
	my(%counts,%events);
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next unless $F[3] eq "nco";
		#next unless $F[3] eq "sco";

		# Skip one of the two discontinous NCO events
		next if $F[1] eq "cs12.16" && $F[2] == 12698579;
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

	print "Chromatids with 2 or more NCOs: $ncoCounts\n";
	print "NCO events < 4Mb from each other: $closeNCOs\n";

#---------------------------------------------------#
# Option 6: Count chromatids with both a CO and NCO #
#---------------------------------------------------#

} elsif ($opts{'a'} == 6) { 
	my(%counts,%events);
	open INF,"./CO_NCO_data/co-ncoDetail.r6.tsv" or die "Can't open file co-ncoDetail.r6.tsv: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;

		# Skip one of the two discontinous NCO events
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
	print "Chromatids with at least one nco and CO event:  $countEventsOnSameChromatid\n";
	print "Chromatids with a NCO event within 4Mb of a CO: $closeNCOs\n";
	print "Average distance between NCO and CO event:      $aveDistBetween\n";
}

# Done
