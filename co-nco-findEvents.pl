#!/usr/bin/perl

use strict;
use Getopt::Std;

my $minCNSBaseQuality	= 60;
my $minVCFScore		= 200; #only use SNPs with scores greater than or equal to this number
my $maxIndelLength	= 1; #ignore indels greater than or equal to this number
my $skipIndels		= 1; # 1 for yes, 0 for no
my $minParentalDepth 	= 20;
my $minChildDepth	= 8;

my $repeatMasker	= "rmsk.txt";
my $chromSizes		= "dm6.chrom.sizes";

my $parentA	= "w1118-parent";
my $parentB	= "CantonS-parent";

my @chrList = qw/ chrX chr2L chr2R chr3L chr3R /; # add chr4 if you'd like to look at that

## Command-line options
my %opts;
getopts('ehls:c:', \%opts); # options as above. Values in %opts

## Usage Output
if ($opts{'h'}) {
        print "
	This script identifies CO and NCO events in child .vcf files from two known
	parental genotypes. Parental files must exist in a parental/ directory,
	chid data must exist in a offspring/ directory.
	
	Required:

		None.

	Optional:

		-e If co-nco-uniqueParentalVariants.tsv already exists you can use it 
		   and not re-make it every time.

		-l Populate or add to the co-nco-PositionsToSkip.tsv file.

		-s Specific stock you want to check.

		-c Chromosome: X, A(utosome), or B(oth). Default B(oth).

		-h This helpful help.
        \n";
        exit 0;
}

## Subroutine to check memory usage - Doesn't work on every OS, notably OSX 
sub memUsage {
        my $gb;
        if (-e "/proc/$$/status") {
                open DATA, "< /proc/$$/status" or die "Unable to read /proc/$$/status: $!\n";
                local $/;
                <DATA> =~ m/^VmSize:\s+(.*?)$/m;
                my($kb) = $1 =~ /(\d+)/;
                $gb = sprintf("%0.1f", $kb / 1000 / 1000);
                $gb .= " Gigabytes";
        } else {
                $gb = "unable to check on this OS\n";
        }
        return $gb;
}

my %qualityScores;
foreach (0..90) {
        $_ += 33;
        my $ascii = chr($_);
        $_ -= 33;
        $qualityScores{$ascii} = $_;
}

## Gather system data
my $pwd = `pwd`;
chomp($pwd);

## Print out all variables for the user
print "[".localtime(time)."] \n";
print "[".localtime(time)."]       Script begin. Variables:\n";
printf "[".localtime(time)."] %30s %8d\n", "minimum cns base quality", $minCNSBaseQuality; 
printf "[".localtime(time)."] %30s %8d\n", "minimum VCF score", $minVCFScore; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for parents", $minParentalDepth; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for child", $minChildDepth; 
if ($skipIndels == 0) {
	printf "[".localtime(time)."] %30s %8d\n", "maximum INDEL length", $maxIndelLength; 
} else {
	printf "[".localtime(time)."] %30s\n", "Skipping all indels"; 
}
print "[".localtime(time)."] \n";

## Grab chromosome names and sizes from the chrom.sizes file
my %chromosomeSizes;
open INF, "$pwd/refGenome/$chromSizes" or die "Can't open $chromSizes: $!";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	# format is: <chromosome>\t<size in nt>
	$chromosomeSizes{$F[0]} = $F[1];
}
close INF;

if (!$opts{'e'}) {  # -e flag is to use existing co-nco-uniqueParentalVariants.tsv file
	## Open repeatmasker file
	my %repeats;
	print "[".localtime(time)."] Getting repeats from $repeatMasker.\n";
	open INF,"$pwd/refGenome/$repeatMasker";
	while (<INF>) {
		next if $_ =~ /Simple_repeat/;
		my(@F) = split /\t/, $_;
		foreach my $id ($F[6]..$F[7]) {
			$repeats{$F[5]}{$id} = 1;
			my $chr = $F[5];
			#$chr =~ s/chr//;
			$repeats{$chr}{$id} = 1;
		}
	}
	close INF;
	print "[".localtime(time)."] Repeats gathered.\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";

	## Grab parental VCF data
	#my $parentCount = 0;
	my %parentalSNPs;
	my %variantPositions;
	print "[".localtime(time)."] Getting parental data.\n";
	foreach my $parent ($parentA, $parentB) {
		my($SNPcount,$skipRepeat,$skipScore,$skipHet,$skipAlt,$skipINDEL,$skipDepth,$indelCount) = (0,0,0,0,0,0,0,0);
		print "[".localtime(time)."] - Getting vcf data from $parent.\n";
		open INF,"$pwd/parental/$parent/$parent.vcf" or die "Can't open $pwd/parental/$parent/$parent.vcf: $!";
		while (<INF>) {
			my(@F) = split /\t/, $_;
			my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);

			next unless $chr =~ /^(chrX|chr2L|chr2R|chr3L|chr3R)$/; # Drosophila specific
			$skipRepeat++ if $repeats{$chr}{$id} == 1;	# skip if in repetative region
			next if $repeats{$chr}{$id} == 1;	
			$skipScore++ if $vcfScore < $minVCFScore;	# skip based on vcf score
			next unless $vcfScore >= $minVCFScore;	
                	$skipHet++ if $F[9] =~ /0\/1/; 			# skip it if it's het when looking at parents
                	next if $F[9] =~ /0\/1/; 	
                	$skipAlt++ if $alt =~ /\,/; 			# skip if the alt allele has two bases
                	next if $alt =~ /\,/; 		

			my $isIndel = 0;
		   	$isIndel = 1 if $_ =~ /INDEL/;
			next if $isIndel == 1 && $skipIndels == 1;

			my $indelLength = 0;
			if ($isIndel == 1) {
				my $refLength = length($ref);
				my $altLength = length($alt);

				$indelLength = abs($refLength - $altLength);
			}
			$skipINDEL++ if $indelLength >= $maxIndelLength;
			next if $indelLength >= $maxIndelLength;

			# calculate depth
			$_ =~ /DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;/;
			my $refDepth = $1 + $2; # this should actually always be 0 for the parents
			my $altDepth = $3 + $4;
			my $totalDepth = $refDepth + $altDepth;
			$skipDepth++ if $totalDepth <= $minParentalDepth;
			next if $totalDepth <= $minParentalDepth;

			#print "$chr\t$id\t$vcfScore\t$refDepth\t$altDepth\t$totalDepth\t$ref\t$alt\n";

			$parentalSNPs{$parent}{$chr}{$id} = "$ref|$alt|$vcfScore|$totalDepth";
			$variantPositions{$chr}{$id} = 1;
			++$SNPcount unless $_ =~ /INDEL/;
			++$indelCount if $_ =~ /INDEL/;
		}
		close INF;

		my $total = $SNPcount + $indelCount;

		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t repeat", $skipRepeat; 
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t score < $minVCFScore", $skipScore;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t het SNP", $skipHet;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t two alt alleles", $skipAlt;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t INDEL >= $maxIndelLength", $skipINDEL;
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t depth <= $minParentalDepth", $skipDepth;
		printf "[".localtime(time)."] %30s %8d\n", "Total INDELs", $indelCount;
		printf "[".localtime(time)."] %30s %8d\n", "Total SNPs", $SNPcount;
		printf "[".localtime(time)."] %30s %8d\n", "Total variants", $total;
		print "[".localtime(time)."] \n";
	}
	print "[".localtime(time)."] Done gathering parental data.\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
	print "[".localtime(time)."] \n";

	## Grab parental consensus sequence data
	my %cns;
	print "[".localtime(time)."] Opening consensus file for each parental line - can be slow.\n";
	foreach my $parent ($parentA, $parentB) {
		print "[".localtime(time)."] - Opening consensus files for $parent.\n";
		foreach my $chr (@chrList) {
			my $cns = "$pwd/parental/$parent/$parent.cns.$chr";
        		open INF,"$cns" or die "Can't open file $cns: $!\n";
			while (<INF>) {
				chomp($_);
				my(@F) = split /\t/, $_;
			
				$F[1] =~ tr/[a-z]/[A-Z]/;
                		$cns{$parent}{$chr}{$F[0]} = "$F[1]|$F[2]";
                	}
			close INF;
        	}
	}
	print "[".localtime(time)."] Done collecting consensus sequence. Memory usage ".memUsage()."\n";
	print "[".localtime(time)."] \n";

	## Identify positions that differ between $parentA and $parentB
	my %skip;
	my %diffVariants;
	my %diffVariantCounts;
	my $output = "Chr\tID\tRef\twcns\tScore\twvcf\tScore\tDepth\tcscns\tScore\tcsvcf\tScore\tDepth\n";
	foreach my $chr (keys %variantPositions) {
		foreach my $id (sort {$a<=>$b} keys %{$variantPositions{$chr}}) {
			my($parentAcns,$parentAcnsScore) = $cns{$parentA}{$chr}{$id} =~ /(\w+)\|(\d+)/;
			my($parentBcns,$parentBcnsScore) = $cns{$parentB}{$chr}{$id} =~ /(\w+)\|(\d+)/;

			$skip{'missingCNS'}++ if !$parentAcns || !$parentBcns;
			next if !$parentAcns || !$parentBcns;

			$skip{'LowCNSScore'}++ if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;
			next if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;

			my($refA,$parentAalt,$AvcfScore,$AvcfDepth) = $parentalSNPs{$parentA}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;
			my($refB,$parentBalt,$BvcfScore,$BvcfDepth) = $parentalSNPs{$parentB}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;

			my $ref = $refA;
		   	$ref = $refB if !$ref;

			$parentAalt 	= "."  if !$parentAalt; 
			$AvcfScore 	= "..." if $parentAalt eq '.'; 
			$AvcfDepth	= ".." if $parentAalt eq '.'; 
			$parentBalt 	= "."  if !$parentBalt; 
			$BvcfScore 	= "..." if $parentBalt eq '.';
			$BvcfDepth 	= ".." if $parentBalt eq '.';

			$skip{'IdenticalAlt'}++ if $parentAalt eq $parentBalt;
			next if $parentAalt eq $parentBalt;

			$diffVariants{$chr}{$id} = "$parentAcns\t$parentAcnsScore\t$parentAalt\t$AvcfScore\t$AvcfDepth\t$parentBcns\t$parentBcnsScore\t$parentBalt\t$BvcfScore\t$BvcfDepth";
			$output .= "$chr\t$id\t$ref\t$diffVariants{$chr}{$id}\n";
			$diffVariantCounts{$chr}++;
	
			# need to add depth of coverage filter at the cns positions
		}
	}

	# deleting structures not needed. See this thread: http://www.perlmonks.org/?node_id=182343
	undef %parentalSNPs;
	undef %cns;

	open OUTF,">$pwd/co-nco-uniqueParentalVariants.tsv";
	print OUTF $output;
	close OUTF;

	undef $output;

	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t missing base", $skip{'missingCNS'}; 
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t low CNS score", $skip{'LowCNSScore'}; 
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t identical variants", $skip{'IdenticalAlt'};

	printf "[".localtime(time)."] \n";
	printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
	my $totalCount;
	foreach my $chr (keys %diffVariantCounts) {
		#my $tmpChr = "chr".$chr;
		my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $diffVariantCounts{$chr});
		printf "[".localtime(time)."] %15s %8d - 1 snp / %4d bp\n", $chr, $diffVariantCounts{$chr}, $snpPerBP; 
		$totalCount += $diffVariantCounts{$chr};
	}
	printf "[".localtime(time)."] %30s %8d\n", "Total Count", $totalCount;
	print "[".localtime(time)."] \n";
	print "[".localtime(time)."] Unique variants saved in co-nco-uniqueParentalVariants.tsv\n";
	print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
	print "[".localtime(time)."] \n";
} else {
	die "Error: co-nco-uniqueParentalVariants.csv doesn't exist! Quitting.\n" if !-e "$pwd/co-nco-uniqueParentalVariants.tsv";
	print "[".localtime(time)."] -e flag passed, co-nco-uniqueParentalVariants.tsv exists.\n";
}

#exit 0;

# open co-nco-uniqueParentalVariants.csv
my(%parental,%SNP,$countVariants,%countSNPs);
print "[".localtime(time)."] Opening co-nco-uniqueParentalVariants.tsv.\n";
open INF,"$pwd/co-nco-uniqueParentalVariants.tsv" or die "Can't open co-nco-uniqueParentalVariants.tsv: $!";
while (<INF>) {
	my(@F) = split /\t/, $_;
	next unless $F[1] =~ /[0-9]/;

	my($chr,$id,$ref,$wcns,$wvcf,$cscns,$csvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

	next if $wvcf ne $wcns && $wvcf !~ /\./;
	next if $csvcf ne $cscns && $csvcf !~ /\./;
	next if $wcns !~ /[A|G|C|T]/ || $cscns !~ /[A|G|C|T]/;
	next if $wcns eq $cscns;

	$SNP{$chr}{$id} = $ref;
	$parental{'w'}{$chr}{$id} = $wcns;
	$parental{'cs'}{$chr}{$id} = $cscns;
	++$countVariants;
	$countSNPs{$chr}++;
}
print "[".localtime(time)."] $countVariants total variants identified.\n";
printf "[".localtime(time)."] \n";
printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
foreach my $chr (@chrList) {
	my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $countSNPs{$chr});
	printf "[".localtime(time)."] %17s %10d - 1 snp / %4d bp\n", $chr, $countSNPs{$chr}, $snpPerBP; 
}
print "[".localtime(time)."] \n";
print "[".localtime(time)."] Memory usage: ".memUsage()."\n";
print "[".localtime(time)."] \n";


######
open INF,"./co-nco-PositionsToSkip.tsv";
while (<INF>) {
	my(@F) = split /\t/, $_;
	my($chr,$count) = ($F[0],$F[2]);
	my($id1,$id2) = $F[1] =~ /(\d+)\-(\d+)/;

	next if $count <= 2;
	
	foreach ($id1..$id2) {
		#print "Skipping: $chr\t$_\n";
		$parental{'w'}{$chr}{$_} = undef;
		$parental{'cs'}{$chr}{$_} = undef;
	}
}
close INF;

my($sep,$stockCount,%eventCount,%eventStocks);

$sep= "+-----+------++-----------+-----+----++-----------+-----+----+-----------+-----++-----------+-----+----++----------+----------++-------------+-------------+";
#print "$sep\n";
#print "| Stk |  Chr ||   Last ID | Bas | Pa ||        ID | Bas | Pa |        ID | Bas ||   Next ID | Bas | Pa ||     fGap |     bGap || Lst->CurSNP | SNPsInEvent | \n";
#print "Chr\tStock\tLastParent\tRange\t\n";


my @files = `ls -1 $pwd/offspring`;
foreach my $stock (@files) {
	#next unless $stock =~ /w[1-9]\./ || $stock =~ /cs[1-9]\./;
	next unless $stock =~ /cs12\.\d/; # || $stock =~ /w13\.\d/;

	#next if $stock =~ /^(w15.1|w15.10|w15.11|w15.12)$/;

	chomp($stock);
	print "[".localtime(time)."] $stock\n";
	$stockCount++;

	my(%vcfData,%allVCFData);

	# open VCF file, store data in %vcf
	next if !-e "$pwd/offspring/$stock/$stock.vcf";
	open INF,"$pwd/offspring/$stock/$stock.vcf" or die "Can't open $pwd/offspring/$stock/$stock.vcf: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);
		$allVCFData{$chr}{$id} = 1;

		next unless $chr =~ /^(chrX|chr2L|chr2R|chr3L|chr3R|chr4)$/; # Drosophila specific
		next unless $parental{'w'}{$chr}{$id} && $parental{'cs'}{$chr}{$id};
		next unless $vcfScore >= $minVCFScore;	
               	next if $alt =~ /\,/; 		

		my $isIndel = 0;
	   	$isIndel = 1 if $_ =~ /INDEL/;
		next if $isIndel == 1 && $skipIndels == 1;

		my $type = "het";
		   $type = "hom" if $F[9] =~ /1\/1/;

               	next if $type eq "het" && $chr eq "chrX"; 	# skip X chr calls that are het

		my $parent = "neither";
		
		if ($type eq "hom" ) {
			$parent = "w" if $alt eq $parental{'w'}{$chr}{$id};
			$parent = "cs" if $alt eq $parental{'cs'}{$chr}{$id};
		} elsif ($type eq "het") {
			$parent = "both" if $alt eq $parental{'w'}{$chr}{$id} && $ref eq $parental{'cs'}{$chr}{$id};
			$parent = "both" if $ref eq $parental{'w'}{$chr}{$id} && $alt eq $parental{'cs'}{$chr}{$id};
		}

		$vcfData{$chr}{$id} = "$parent|$type|$ref|$alt";
		#print "$vcfData{$chr}{$id}\n" if $id < 1000000 && $chr eq "chr2L";
	}

	foreach my $chr (@chrList) {
		#next if $chr eq "chrX";
		foreach my $id (keys %{$SNP{$chr}}) {
			next if $vcfData{$chr}{$id};
			next if $allVCFData{$chr}{$id} == 1;

			my $parent = "neither";
			$parent = "w"  if $SNP{$chr}{$id} eq $parental{'w'}{$chr}{$id};
			$parent = "cs" if $SNP{$chr}{$id} eq $parental{'cs'}{$chr}{$id};

			$vcfData{$chr}{$id} = "$parent|ref|$SNP{$chr}{$id}|N";
		}


		#print "[".localtime(time)."] $stock | $chr \n";
		my($lastParent,$lastParentChange,$lastID,$snpCount);
		my($changeList);
		foreach my $id (sort {$a<=>$b} keys %{$vcfData{$chr}}) {
			my($parent,$type,$ref,$alt) = split /\|/, $vcfData{$chr}{$id};
			#print "$vcfData{$chr}{$id}\n" ; #if $id == 7023293;

			next if $parent eq "neither";

			if (!$lastParent) {
				$lastParent = $parent;
				$lastParentChange = $id;
				$snpCount = 1;
				$lastID = $id;
			} elsif ($parent eq $lastParent) {
				$snpCount++;
				$lastID = $id;
			} elsif ($parent ne $lastParent) {
				my($depth,$allele,$alleles) = getAllele($stock,$chr,$id);

				my $tmpParent = "neither";
				   $tmpParent = "w" if $allele eq $parental{'w'}{$chr}{$id} && $alleles !~ /\|\w\|/;
				   $tmpParent = "cs" if $allele eq $parental{'cs'}{$chr}{$id} && $alleles !~ /\|\w\|/;

				if ($alleles =~ /(\w)\|(\w)\|/) {
					my($a,$b) = ($1,$2);
					$tmpParent = "both" if $a eq $parental{'w'}{$chr}{$id} && $b eq $parental{'cs'}{$chr}{$id};
					$tmpParent = "both" if $b eq $parental{'w'}{$chr}{$id} && $a eq $parental{'cs'}{$chr}{$id};
				}

				#print "$stock\t$chr\t$id\t$depth\t$allele\t$tmpParent\t$parent\n" if $id == 24865481;
				if ($depth > $minChildDepth && $tmpParent eq $parent) {
					if ($changeList) {
						my $gap = $id - $lastID;
						my $lastGap = $lastID - $lastParentChange;
						print "$changeList-$lastID\t$lastGap\t$snpCount\t-> $parent $id $gap\n";
						$changeList = undef;

						if ($snpCount < 100) {
							$eventCount{$chr}{"$lastParentChange-$lastID"}++;
							$eventStocks{$chr}{"$lastParentChange-$lastID"} .= "$stock,";
						}
						
					} else {
						my $gap = $id - $lastID;
						my $lastGap = $lastID - $lastParentChange;
						$changeList = "$chr\t$stock\t$lastParent\t$lastParentChange-$lastID\t$snpCount\t$lastGap\t->\t$gap\t$parent\t$chr:$id";
					}

					$lastID = $id;
					$lastParent = $parent;
					$lastParentChange = $id;
					$snpCount = 1;
				}
			}
		}
		if ($changeList =~ /$chr/) {
			my $lastGap = $lastID - $lastParentChange;
			print "$changeList-$lastID\t$lastGap\t$snpCount\t-> END_OF_CHR\n";
		}
	}
}

if ($opts{'l'}) {
	open OUTF,">cogcPositionsToSkip.tsv";
	print OUTF "Chr\tRange\tCount\tStocks Seen In (stocks checked: $stockCount)\n";
	foreach my $chr (keys %eventCount) {
		foreach my $id (sort {$a<=>$b} keys %{$eventCount{$chr}}) {
			my $count = $eventCount{$chr}{$id};
			#next unless $count > 3;
			my($id1,$id2) = $id =~ /(\d+)\-(\d+)/;
			my $gap = $id2 - $id1;
			print OUTF "$chr\t$id\t$count\t$gap\t$eventStocks{$chr}{$id}\n";
		}
	}
	close OUTF;
}

#####################

sub getAllele {
	my($stock,$chr,$id) = ($_[0],$_[1],$_[2]);

	my %alleleCount;

	foreach my $chunk (split /\n/, `samtools view $pwd/offspring/$stock/$stock.bam $chr:$id-$id`) {
		my(@F) = split /\t/, $chunk;
		next if $F[4] < 60;
		my($splitInfo,$originalSequence,$originalQuality) = ($F[5],$F[9],$F[10]);
	
       		my @seq;
       		foreach (split //, $originalSequence) {
       			push(@seq,$_);
       		}
       		$originalSequence = undef;
       		my $count = 0;
       		foreach (split //, $originalQuality) {
       			my $score = $qualityScores{$_};
       			my $base = $seq[$count];
       			$base = "n" if $qualityScores{$_} < 15;
       			$originalSequence .= $base;
       			++$count;
       		}

        	my $unknownAction = 0;
		my $seq;
        	foreach (split /([0-9]+[A-Z])/, $splitInfo) {
        		my($len,$action) = $_ =~ /([0-9]+)([A-Z])/;
        		next unless $len =~ /[0-9]/;

        		if ($action =~ /M/) {
        			$originalSequence =~ /^(.{$len})/;
        			$seq .= $1;
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /D/) {
        			foreach (1..$len) {
        				$seq .= ".";
        			}
        		} elsif ($action =~ /I/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /S/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} else {
				$unknownAction = 1;
			}
		}

		next if !$seq;
		my $length = length($seq);
		my $startID = $F[3];
		my $endID = $F[3] + $length - 1;

		$startID = $endID if $startID > $id;

		my $currID = $startID;
		foreach (split //, $seq) {
			#print "$currID\t$_\n" if $currID == 7023293;
			$alleleCount{$_}++ if $currID == $id;
			$currID++;
		}

	}

	my $readCount = $alleleCount{'A'} + $alleleCount{'G'} + $alleleCount{'C'} + $alleleCount{'T'};
	my($allele,$alleles);
	if ($readCount > 0) {
		my $tmpReadCount = $readCount / 2;
		$allele = "A" if $alleleCount{'A'} >= $tmpReadCount;
		$allele = "G" if $alleleCount{'G'} >= $tmpReadCount;
		$allele = "C" if $alleleCount{'C'} >= $tmpReadCount;
		$allele = "T" if $alleleCount{'T'} >= $tmpReadCount;

		my $alleleFreq = sprintf("%0.0f", (($alleleCount{$allele} / $readCount ) * 100));

		if ($alleleFreq >= 30 && $alleleFreq <= 70) {
			$alleles .= "A|" if $alleleCount{'A'} >= 2;
			$alleles .= "G|" if $alleleCount{'G'} >= 2;
			$alleles .= "C|" if $alleleCount{'C'} >= 2;
			$alleles .= "T|" if $alleleCount{'T'} >= 2;
		} elsif ($alleleFreq >= 5 && $alleleFreq <= 95) {
			$allele = "N";
		}

		#print "$alleleFreq\t$allele\t$alleles\t$alleleCount{'A'}, $alleleCount{'G'}, $alleleCount{'C'}, $alleleCount{'T'}\n" if $id == 24865481;
	}

	return($readCount,$allele,$alleles);
}
