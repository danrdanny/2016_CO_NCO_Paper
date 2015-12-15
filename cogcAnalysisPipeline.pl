#!/usr/bin/perl

# Description
# 
# This pipeline takes fastq data for the CO/NCO paper and aligns/analyzes it.
#
# Please use -h to see a full description. 
# Stocks to align/analyze are listed in sampleSheet.tsv. 
# Lines beginning with '#' are ignored.
#
# Written by Danny Miller danrdanny (shift-2) gmail (dot) com

use strict;
use Getopt::Std;

## Global Variables
my $samQuality 		= 30; # passed to bwa, any removes any read with a overall quality score < this number
my $minReadPairs 	= 10; # breakdancer: minimum number of read pairs needed to call a SV

## Paths for all necessary programs
# Realignment steps
my $GATK 		= "/n/local/stage/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar";
my $picard 		= "/n/local/stage/picard/picard_1.129/picard.jar";
#Identify split/discordant reads
my $samblaster 		= "/home/dem/bin/samblaster/samblaster";
#SV detection with Breakdancer
my $breakdancercfg 	= "/home/dem/bin/breakdancer/perl/bam2cfg.pl";
my $breakdancer 	= "/home/dem/bin/breakdancer/bin/breakdancer-max";
#SNP Calling
my $samtools 		= "/n/local/bin/samtools-0.1.19";
my $bcftools 		= "/n/local/bin/bcftools-0.1.19";
my $vcfutils 		= "/n/local/stage/samtools/samtools-0.1.19/bcftools/vcfutils.pl";
my $vcfMerge		= "/n/local/bin/vcf-merge";
my $vcftools		= "/n/local/bin/vcftools";
#Java variables
my $javaFlags		= "-Xmx4g";

#Other programs used, should generally be in $PATH
my @programs = qw/ rm mv cat mkdir wget zcat cp echo gzip tar awk bgzip tabix java twoBitToFa bwa genomeCoverageBed bedtools bedGraphToBigWig /;

foreach my $program (@programs) {
	my $progPath = executeCommand("which $program 2>&1");
	chomp($progPath);

	if ($progPath =~ /(no $program|command not found)/ig) {
		print "Path to $program not found. Please add to your \$PATH. Exiting.\n";
		exit 0;
	}
}

#------# Nothing should need to be changed below this line to run the program successfully #------#

## Global Variables
my $repeatMasker	= "rmsk.txt";
my $repeatMaskerGzName	= "rmsk.txt.gz";
my $repeatMaskerURL	= "http://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz";
my $dm6Ref		= "dm6.fa";
my $ucscRefName		= "dm6.fa.gz";
my $ucscdm6 		= "http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz";
my $chromSizes		= "dm6.chrom.sizes";
my $ucscChromSizes 	= "http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes";
my $genomeFile		= "genome.dm6.txt"; # from github
#my $md5File		= "md5values.individualFiles.dat"; # from github
my $sampleSheet		= "sampleSheet.tsv"; # from github

## Logfile name
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1; 
$mon = "0$mon" if $mon =~ /^\d$/; # this is because I'm OCD and can't stand when the log file names don't look nice
$mday = "0$mday" if $mday =~ /^\d$/;
$hour = "0$hour" if $hour =~ /^\d$/;
$min = "0$min" if $min =~ /^\d$/;
$sec = "0$sec" if $sec =~ /^\d$/;
my $logfile = "fm7.alignment.$year-$mon-$mday\_$hour\:$min\:$sec\.log";

## Command-line options
my %opts;
getopts('t:ah', \%opts); # Values in %opts

## Usage Statement
if ($opts{'h'} || !$opts{'t'}) {
	print "

	This script takes fastq files and aligns them to the dm6 reference genome. This script
	was used to align the data used in the following publication:

	A whole-genome analysis of individual meiotic events in Drosophila melanogaster reveals 
	that interference and the centromere effect are properties unique to each chromosome arm

	All commands are logged in ./log/<timestamp>.log

	Generally, the steps are:
	  - Check for the presence of required files: dm6.fa, dm6.chrom.sizes, rmsk.txt, and
	    genome.dm3.txt. Any file absent will be downloaded automatically.

	  - Check for the presence of sampleSheet.tsv. Any line that begins with a # in this
	    file will be ignored. Format of this file is as follows:
	    <Sample_Name>\\t<Barcode>\\t<LIMS_Order>\\t<Lanes>\\t<Class>

	  - Check for the presence of raw data files. The script will only download files 
	    needed to complete the alignment of stocks listed in sampleSheet.tsv.

	  - Perform alignment with bwa on each stock in sampleSheet.tsv not prefixed with '#'.
	    If .bam file exists in the <Sample_Name> directory this step is skipped.
	    - After merging/sorting read groups are replaced with picard 
	      AddOrReplaceReadGroups.jar
	    - Duplicates are marked with picard MarkDuplicates.jar
	    - Read groups are again replaced with picard AddOrReplaceReadGroups.jar
	    - During the alignment process split and discordant reads are identified and saved
	      in a sam file.

	  - Calculate genome coverage. This step is skipped if <Sample_Name>/<Sample_Name>.bw 
	    exists.

	  - Calculate insert sizes. This step is skipped if 
            <Sample_Name>/<Sample_Name>.insertMetrics.pdf exists.

	  - Calculate alignment statistics. This step is skipped if 
	    <Sample_Name>/<Sample_Name>.alignStats exists.

	  - Call SNPs with samtools. This step is skipped if <Sample_Name>/<Sample_Name>.vcf
	    exists.

	  - Collect alignment statistics. This step is skipped if 
	    <Sample_Name>/<Sample_Name>.alignmentStatistics.tsv exists.

	  - Additional steps are outlined in the optional flags below.

	Required flags:
	-t  Threads you want any multi-threaded program to run with.

	Optional flags:
	-a  Make common alignmentStatistics.tsv file.
            Specifically, this will collect data from each <Sample_Name>/alignmentStatistics.tsv 
	    file and pool it into an out_alignmentStatistics.tsv file at the root directory. It 
	    will operate over stocks listed in sampleSheet.tsv.

	-h  This helpful help.
	\n";
	exit 0;
}

## Get our current working directory
my $pwd = `pwd`; chomp($pwd);
my $refLoc = $pwd."/refGenome";

## Threads to run multithreadded programs with
my $threads = $opts{'t'};

## Two subroutines are used, one to run commands, the second to log all activities
sub executeCommand {
	open LOGF, ">>$pwd/out_log/$logfile";
	print LOGF "[".localtime()."] CMD: $_[0]\n";
	close LOGF;
	my $output = `$_[0]`;
	return($output);
}

sub logData {
	print "[".localtime()."] $_[0]\n" if $_[1] eq "print";
	open LOGF, ">>$pwd/out_log/$logfile";
	print LOGF "[".localtime()."] LOG: $_[0]\n";
	close LOGF;
	return 1;
}
## End Subroutines 

## Make our log directory
executeCommand("mkdir -p out_log");

# -----------------------------------------------------------------------#
# Check if required files exist. 
# Some are downloaded if missing, others should be pulled in from github. 
# There are 6 required files.
# -----------------------------------------------------------------------#

# 1. dm6.fa - reference genome
if (!-e "$refLoc/$dm6Ref") {
	print "Downloading reference genome from $ucscdm6\n";
	executeCommand("wget $ucscdm6");
	executeCommand("gzip -d $ucscRefName");

	executeCommand("bwa index $dm6Ref");
	executeCommand("java -jar $picard/CreateSequenceDictionary.jar R=$dm6Ref O=dm6.dict");
	executeCommand("$samtools faidx $dm6Ref");
	executeCommand("mv dm6.* $refLoc");
}
die "Could not locate reference genome $dm6Ref" if !-e "$refLoc/$dm6Ref";
logData("Identified reference genome $dm6Ref",'print');

# 2. chrom.sizes - to calculate coverage per arm
if (!-e "$refLoc/$chromSizes") {
	print "Downloading chromSizes files from $ucscChromSizes\n";
	executeCommand("wget $ucscChromSizes");
	executeCommand("mv $ucscChromSizes $refLoc");
}
die "Could not locate chrom.sizes file" if !-e "$refLoc/$chromSizes";
logData("Identified chrom.sizes file $chromSizes",'print');

# 3. genome.dm6.txt - for bedtools to calculate coverage
die "Could not locate $genomeFile file" if !-e "$refLoc/$genomeFile";
logData("Identified file $genomeFile",'print');

# 4. sampleSheet.tsv
die "Could not locate $sampleSheet file" if !-e $sampleSheet;
logData("Identified $sampleSheet file",'print');

# 5. md5values of individual files
#die "Could not locate $md5File file" if !-e $md5File;
#logData("Identified md5 summary file $md5File");

# 6. repeat masker file - to ignore repetitive regions
if (!-e "$refLoc/$repeatMasker") {
	print "Downloading repeat masker file from $repeatMaskerURL\n";
	executeCommand("wget $repeatMaskerURL");
	executeCommand("gzip -d $repeatMaskerGzName");
	executeCommand("mv $repeatMasker $refLoc");
}
die "Could not locate repeat masker file" if !-e "$refLoc/$repeatMasker";
logData("Identified repeat masker file $repeatMasker",'print');


## Grab chromosome names and sizes from the chrom.sizes file
my %chromosomeSizes;
open INF, "$refLoc/$chromSizes" or die "Can't open $chromSizes: $!";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	# format is: <chromosome>\t<size in nt>
	$chromosomeSizes{$F[0]} = $F[1];
}
close INF;

# ---------------------------------------------------------------#
# Check if all necessary fastq data is present, if not, download #
# ---------------------------------------------------------------#

## open md5value sheet to get fastq files associated with each stock
#my %fastqFiles;
#open INF,"$pwd/$md5File" or die "Can't open $md5File: $!";
#while (<INF>) {
#	chomp($_);
#	my($md5,$lims,$file) = $_ =~ /(\w+)\s+(MOLNG-\d+)\/(\w+)\/(\S+)/;
#	$fastqFiles{$lims}{$file} = $md5;
#}
#close INF;

## Open sampleSheet.tsv to get information about stocks/barcodes
my %stocksToAlign;
#my %limsOrdersNeeded;
my $numStocksToAlign;
#my $rawDataMissing = 0;
open INF,"$pwd/$sampleSheet" or die "Can't open $sampleSheet: $!";
while (<INF>) {
	chomp($_);
	next if $_ =~ /^#/; # skip any line that starts with #
	next unless $_ =~ /[a-z0-9]/i;
	logData("Stock: $_");
	my(@F) = split /\t/, $_;

	my $sampleName 	= $F[0];
	my $barcode	= $F[1];
	my $LIMS_Order	= $F[2];
	my($laneStart,$laneEnd) = $F[3] =~ /(\d+)\-(\d+)/;

	print "[".localtime(time)."] $sampleName | Checking for files in $LIMS_Order ... ";
	foreach my $lane ($laneStart..$laneEnd) {
		print "$LIMS_Order/s_$lane\_1_$barcode.fastq.gz missing - exiting\n" if !-e "$pwd/$LIMS_Order/s_$lane\_1_$barcode.fastq.gz";
		die if !-e "$pwd/$LIMS_Order/s_$lane\_1_$barcode.fastq.gz";

		print "$LIMS_Order/s_$lane\_1_$barcode.fastq.gz missing - exiting\n" if !-e "$pwd/$LIMS_Order/s_$lane\_2_$barcode.fastq.gz";
		die if !-e "$pwd/$LIMS_Order/s_$lane\_2_$barcode.fastq.gz";

		#if (!-e "$pwd/$LIMS_Order/$file") {
		#	print "\n";
		#	logData("Downloading data for $LIMS_Order/$flowcell/$file",'print');
		#	executeCommand("wget $ftpHost/$LIMS_Order/$flowcell/$file");
		#	executeCommand("$mkdir -p $LIMS_Order/$flowcell/");
		#	executeCommand("$mv $file $LIMS_Order/$flowcell/");
		#}
	}
	print " OK\n";

	$stocksToAlign{$sampleName} = $_;
	++$numStocksToAlign;
}
close INF;

logData("All necessary data files present",'print');

#-----------------#
# Begin Alignment #
#-----------------#

## Loop over the stocks we should be aligning. 
foreach my $newStockName (keys %stocksToAlign) {
	print "[".localtime(time)."] $newStockName | ----- Start -----\n";

	my(@F) = split /\t/, $stocksToAlign{$newStockName};

	my $targetDir = "$pwd/offspring/$newStockName";
	executeCommand("mkdir -p $targetDir");

	my $barcode 	  	= $F[1];
	my $LIMS_Order	  	= $F[2];
	my($laneStart,$laneEnd) = $F[3] =~ /(\d+)\-(\d+)/;
	my $stockType		= $F[4];

	my $bamFileList; # list of bam files so they can be merged, then deleted
	my $finalBam = $newStockName . ".bam";
	my $splitDiscFinal = "$newStockName.split.discordant.sam";
	my $unmappedFinal = "$newStockName.unmapped.sam";
	my $fastqLocation = "$pwd/$LIMS_Order";

	# does the fastq directory exist? If not, die.
	die "fastq files not found at $fastqLocation\n" if !-e $fastqLocation;

 	# First check to see if a final .bam file exists, if it does, we can skip this
 	if (!-e "$targetDir/$finalBam") {
		my $fastqCount = 0; # 
		foreach my $lane ($laneStart..$laneEnd) {
			my $pair1 = "s_".$lane."_1_$barcode.fastq.gz";
			my $pair2 = "s_".$lane."_2_$barcode.fastq.gz";

			$fastqCount++;
			logData("Running analysis on pair1: $pair1, pair2: $pair2");

			# Align with BWA
			print "[".localtime(time)."] $newStockName | RUN: Aligning $pair1 and $pair2 with bwa\n";
			executeCommand("bwa mem -t $threads $refLoc/$dm6Ref $fastqLocation/$pair1 $fastqLocation/$pair2 > $targetDir/$barcode.$fastqCount.sam 2>/dev/null");

			# convert sam to bam
			executeCommand("$samtools view -q $samQuality -bS $targetDir/$barcode.$fastqCount.sam > $targetDir/$barcode.$fastqCount.bam 2>>$pwd/out_log/stderr");
			$bamFileList .= "$targetDir/$barcode.$fastqCount.bam ";
	
			# run samblaster
			print "[".localtime(time)."] $newStockName | RUN: running samblaster on $barcode.$fastqCount.sam\n";
			my($discordant,$split,$unmapped) = ("$targetDir/discordant.tmp","$targetDir/split.tmp","$targetDir/unmapped.tmp");
			executeCommand("$samblaster -i $targetDir/$barcode.$fastqCount.sam -d $discordant -s $split -u $unmapped -o /dev/null 2>>$pwd/out_log/stderr");
			executeCommand("cat $discordant | grep -v \"^\@\" >> $targetDir/$splitDiscFinal");
			executeCommand("cat $split | grep -v \"^\@\" >> $targetDir/$splitDiscFinal");
			executeCommand("cat $unmapped >> $targetDir/$unmappedFinal");
			executeCommand("rm -f $discordant $split $unmapped");

			# remove sam files
			executeCommand("rm -f $targetDir/$barcode.$fastqCount.sam");
		}

		executeCommand("gzip $targetDir/$splitDiscFinal");
		executeCommand("rm -f $targetDir/$unmappedFinal");

		print "[".localtime(time)."] $newStockName | RUN: Merging and sorting bam files.\n";
		executeCommand("$samtools merge -f $targetDir/$barcode.merge.bam $bamFileList");
		executeCommand("rm -f $bamFileList");
		executeCommand("$samtools sort -\@ $threads $targetDir/$barcode.merge.bam $targetDir/$barcode.sorted 2>>$pwd/out_log/stderr");

		# Picard step 1, AddOrReplaceReadGroups
		print "[".localtime(time)."] $newStockName | RUN: Picard step 1: AddOrReplaceReadGroups\n";
		executeCommand("java -jar $picard AddOrReplaceReadGroups I=$targetDir/$barcode.sorted.bam O=$targetDir/$barcode.picardIntermediate.1.bam RGLB=$newStockName RGPL=illumina RGPU=$barcode RGSM=$newStockName CREATE_INDEX=TRUE 2>>$pwd/out_log/stderr");

		# Picard step 2, MarkDuplicates
		print "[".localtime(time)."] $newStockName | RUN: Picard step 2: MarkDuplicates\n";
		executeCommand("java -jar $picard MarkDuplicates I=$targetDir/$barcode.picardIntermediate.1.bam O=$targetDir/$barcode.picardIntermediate.2.bam METRICS_FILE=$targetDir/$barcode.picardIntermediate.metrics.dat REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE 2>>$pwd/out_log/stderr");

		# GATK step 1, RealignerTargetCreator
		print "[".localtime(time)."] $newStockName | RUN: GATK step 1: RealignerTargetCreator\n";
		executeCommand("java $javaFlags -jar $GATK -T RealignerTargetCreator -I $targetDir/$barcode.picardIntermediate.2.bam -R $refLoc/$dm6Ref -o $targetDir/$barcode.GATKIntermediate.intervals 2>>$pwd/out_log/stderr");

		# GATK step 2, IndelRealigner 
		print "[".localtime(time)."] $newStockName | RUN: GATK step 2: IndelRealigner\n";
		executeCommand("java $javaFlags -jar $GATK -T IndelRealigner -I $targetDir/$barcode.picardIntermediate.2.bam -R $refLoc/$dm6Ref -o $targetDir/$barcode.GATKIntermediate.1.bam -targetIntervals $targetDir/$barcode.GATKIntermediate.intervals 2>>$pwd/out_log/stderr");

		executeCommand("mv $targetDir/$barcode.GATKIntermediate.1.bam $targetDir/$finalBam");
		executeCommand("$samtools index $targetDir/$finalBam");

		# Remove files used/generated during this step
		executeCommand("rm -f $targetDir/$barcode.merge.bam");
		executeCommand("rm -f $targetDir/$barcode.sorted.bam");
		executeCommand("rm -f $targetDir/*picardIntermediate*");
		executeCommand("rm -f $targetDir/*GATKIntermediate*");
		executeCommand("rm -f $targetDir/*.fastq");
	} else {
		print "[".localtime(time)."] $newStockName | $finalBam exists, skipping alignment.\n"; 
	}

	# Create consensus sequence file for parental lines
	if (!-e "$targetDir/$newStockName.cns.fq.gz" && $stockType eq "parent") {
		print "[".localtime(time)."] Creating consensus file for $newStockName.\n";
		executeCommand("$samtools mpileup -uf $refLoc/$dm6Ref $targetDir/$newStockName.bam 2>>$pwd/out_log/stderr | $bcftools view -cg - 2>>$pwd/out_log/stderr | $vcfutils vcf2fq > $targetDir/$newStockName.cns.fq 2>>$pwd/out_log/stderr");

		my $cns = $newStockName.".cns.fq";
		print "[".localtime(time)."] Opening $newStockName consensus file.\n";
        	open INF,"$targetDir/$cns" or die "Can't open file $cns: $!\n";
        	my $chr;
        	my %cnsData;

        	while (<INF>) {
                	chomp($_);

                	if ($_ =~ /^\@(\w+)$/) {
                        	$chr = $1;
                	} elsif ($_ =~ /^(\+)$/) {
                        	$chr .= "|+";
                	} else {
                        	$cnsData{$chr} .= $_;
                	}
        	}

        	foreach my $chr (keys %cnsData) {
                	next if $chr =~ /\|\+/;
			next unless $chr =~ /^(chrX|chr2L|chr2R|chr3L|chr3R|chr4)$/;

                	my %seqData;
                	my $seqDataCount = 0;
                	foreach my $foo (split //, $cnsData{$chr}) {
                        	next unless $foo =~ /\w/;
                        	++$seqDataCount;
                        	$seqData{$seqDataCount} = $foo;
                	}
                	delete($cnsData{$chr});

                	my %fqData;
                	my $fqDataCount = 0;
                	foreach my $foo (split //, $cnsData{"$chr|+"}) {
                        	next unless $foo =~ /\S/;
                        	++$fqDataCount;
                        	my $qual = ord($foo) - 33;
                        	$fqData{$fqDataCount} = $qual;
                	}
                	delete($cnsData{"$chr|+"});

			my $output;
                	foreach my $id (1..$seqDataCount) {
                		#$cns{$parent}{$chr}{$id} = "$seqData{$id}|$fqData{$id}";
                		$output .= "$id\t$seqData{$id}\t$fqData{$id}\n";
                	}
                	
                	open OUTF,">$targetDir/$newStockName.cns.$chr";
                	print OUTF $output;
                	close OUTF;

			executeCommand("gzip $targetDir/$newStockName.cns.$chr");
                }
                %cnsData = undef;
		print "[".localtime(time)."] Done collecting consensus sequence.\n";
		print "[".localtime(time)."] \n";

		executeCommand("gzip $targetDir/$newStockName.cns.fq");
	}

	# Run Breakdancer
	if (!-e "$targetDir/$newStockName.ctx") {
		print "[".localtime(time)."] $newStockName | RUN: Running breakdancer on $newStockName\n";
 		executeCommand("$breakdancercfg -m $targetDir/$finalBam > $targetDir/$newStockName.cfg 2>>$pwd/out_log/stderr");

		# add -l if you are using mate-pair data
 		executeCommand("$breakdancer -r $minReadPairs -d $targetDir/$finalBam $targetDir/$newStockName.cfg > $targetDir/$newStockName.ctx 2>>$pwd/out_log/stderr");
		executeCommand("rm -f $targetDir/$newStockName.cfg");
 	} else {
		print "[".localtime(time)."] $newStockName | OK:  Breakdancer file $newStockName.ctx exists, skipping this step.\n";
	}

	# Calculate genome coverage
	if (!-e "$targetDir/$newStockName.bw") {
		print "[".localtime(time)."] $newStockName | RUN: Calculating genome coverage for $newStockName\n";
		executeCommand("genomeCoverageBed -bga -trackline -trackopts 'name=\"$finalBam\" visibility=2' -ibam $targetDir/$finalBam -g $refLoc/$genomeFile > $targetDir/$newStockName.coverage");
 		executeCommand("bedGraphToBigWig $targetDir/$newStockName.coverage $refLoc/$genomeFile $targetDir/$newStockName.bw");
		executeCommand("rm -f $targetDir/$newStockName.coverage");
	} else {
		print "[".localtime(time)."] $newStockName | OK:  Genome coverage file exists for $newStockName\n";
	}

	# Calc insert size
	if (!-e "$targetDir/$newStockName.insertMetrics.pdf") {
		print "[".localtime(time)."] $newStockName | RUN: Calculating insert size for $newStockName\n";
		executeCommand("java -jar $picard CollectInsertSizeMetrics I=$targetDir/$finalBam O=$targetDir/$newStockName.insertMetrics H=$targetDir/$newStockName.insertMetrics.pdf 2>>$pwd/out_log/stderr");
	} else {
		print "[".localtime(time)."] $newStockName | OK:  Insert metrics file exists for $newStockName\n";
	}

	# Calc Alignment Stats including ave read size
	if (!-e "$targetDir/$newStockName.alignStats") {
		print "[".localtime(time)."] $newStockName | RUN: Calculating alignment statistics for $newStockName\n";
		executeCommand("java -jar $picard CollectAlignmentSummaryMetrics I=$targetDir/$finalBam O=$targetDir/$newStockName.alignStats 2>>$pwd/out_log/stderr");
	} else {
		print "[".localtime(time)."] $newStockName | OK:  Alignment statistics file exists for $newStockName\n";
	}

	# Call snps using samtools
	if (!-e "$targetDir/$newStockName.vcf") {
		print "[".localtime(time)."] $newStockName | RUN: Calling SNPs for $newStockName\n";
		executeCommand("$samtools mpileup -u -f $refLoc/$dm6Ref $targetDir/$finalBam 2>>$pwd/out_log/stderr | $bcftools view -bvcg - > $targetDir/$newStockName.bcf 2>>$pwd/out_log/stderr");
		executeCommand("$bcftools view $targetDir/$newStockName.bcf | $vcfutils varFilter -D500 > $targetDir/$newStockName.vcf");
		executeCommand("rm -f $targetDir/$newStockName.bcf");
	} else {
		print "[".localtime(time)."] $newStockName | OK:  VCF file exists for $newStockName\n";
	}

	#------# Calculate alignment metrics #------#
	if (!-e "$targetDir/$newStockName.alignmentStatistics.tsv") {
		my %alignmentStatistics;
		my $finalOutput;
		# Count # of reads aligned
		open INF,"$targetDir/$newStockName.alignStats";
		while (<INF>) {
			my @F = split /\t/, $_;
			next unless $F[1] =~ /[1-9]/;

			$alignmentStatistics{$newStockName}{'Picard_Read_Length'} = $F[15];
			$alignmentStatistics{$newStockName}{'Picard_First_Of_Pair'} = $F[1] if $F[0] =~ /FIRST_OF_PAIR/;
			$alignmentStatistics{$newStockName}{'Picard_Second_Of_Pair'} = $F[1] if $F[0] =~ /SECOND_OF_PAIR/;
		}
		close INF;

		# Count reads that align to X, 2L, 2R, 3L, 3R, or the 4th 
		print "[".localtime(time)."] $newStockName | STAT: Gathering per-chromosome alignment information\n";
		foreach my $chr ('chrX','chr2L','chr2R','chr3L','chr3R','chr4') {
			my $genomeAligned = executeCommand("$samtools view $targetDir/$newStockName.bam $chr | wc -l");
			chomp($genomeAligned);
			$alignmentStatistics{$newStockName}{$chr} = $genomeAligned;
			#print "[".localtime(time)."]   $chr\t$genomeAligned\n";
		}

		# Count fastq reads
		my $fastqFileCount;
		print "[".localtime(time)."] $newStockName | STAT: Counting total fastq reads\n";
		my $fastqLocation = "$pwd/$LIMS_Order";
		foreach my $lane ($laneStart..$laneEnd) {
			my $fileName = $fastqLocation."/s_".$lane."_1_".$barcode.".fastq.gz";
			my $fastqReadCount = executeCommand("zcat $fileName | echo \$((`wc -l`/4))");
			$fastqReadCount *= 2;
			$alignmentStatistics{$newStockName}{'Total_fastq_reads'} += $fastqReadCount;
			$fastqFileCount += 2;
		}

		# get Insert size
		open INF,"$pwd/offspring/$newStockName/$newStockName.insertMetrics";
		while (<INF>) {
			my(@G) = split /\t/, $_;
			next unless $G[4] > 0;
			$alignmentStatistics{$newStockName}{'Insert_Size'} = sprintf("%0.0f", $G[4]);
		}
		close INF;

		# Output
		$finalOutput .= "$newStockName\t$barcode\t$LIMS_Order\t$laneStart-$laneEnd\t";
		$finalOutput .= "$alignmentStatistics{$newStockName}{'Picard_Read_Length'}\t"; 
		$finalOutput .= "$alignmentStatistics{$newStockName}{'Insert_Size'}\t";
		$finalOutput .= "$fastqFileCount\t";
		$finalOutput .= "$alignmentStatistics{$newStockName}{'Total_fastq_reads'}\t";
		$finalOutput .= "$alignmentStatistics{$newStockName}{'Picard_First_Of_Pair'}\t";
		$finalOutput .= "$alignmentStatistics{$newStockName}{'Picard_Second_Of_Pair'}\t";

		foreach my $chr ('chrX','chr2L','chr2R','chr3L','chr3R','chr4') {
			$finalOutput .= "$alignmentStatistics{$newStockName}{$chr}\t";
		}

		foreach my $chr ('chrX','chr2L','chr2R','chr3L','chr3R','chr4') {
			my $depthOfCoverage = sprintf("%0.0f", ( ($alignmentStatistics{$newStockName}{$chr} * $alignmentStatistics{$newStockName}{'Picard_Read_Length'}) / $chromosomeSizes{$chr}));
			$finalOutput .= $depthOfCoverage . "x\t";
		}

		my $x_vs_2L = sprintf("%0.2f", $alignmentStatistics{$newStockName}{'chrX'} / $alignmentStatistics{$newStockName}{'chr2L'});
		my $fourth_vs_2L = sprintf("%0.2f", $alignmentStatistics{$newStockName}{'chr4'} / $alignmentStatistics{$newStockName}{'chr2L'});

		$finalOutput .= $x_vs_2L . "\t" . $fourth_vs_2L;
		$finalOutput .= "\n";

		# Print summary statistics to a file
		my $finalOutputHeader  = "NewStock\tBarcode\tLIMS_Order\tFlowcell\tReadLength\tAveInsertSize\tFASTQFiles\tTotalFastqReads\t";
   		$finalOutputHeader .= "AlignedFirstOfPair\tAlignedSecondOfPair\t";
   		$finalOutputHeader .= "chrX\tchr2L\tchr2R\tchr3L\tchr3R\tchr4\t";
   		$finalOutputHeader .= "chrXDepth\tchr2LDepth\tchr2RDepth\tchr3LDepth\tchr3RDepth\tchr4Depth\tX:2L\t4:2L";
   		$finalOutputHeader .= "\n";

		open OUTF,">$targetDir/$newStockName.alignmentStatistics.tsv";
		print OUTF $finalOutputHeader;
		print OUTF $finalOutput;
		close OUTF;
	}
	
	print "[".localtime(time)."] $newStockName | ----- Done -----\n";
	print "[".localtime(time)."]\n";
}

# Generate a common alignment statistics file
#if ((!-e "$pwd/alignmentStatistics.tsv" && $opts{'a'}) || $opts{'A'}) {
if ($opts{'a'}) {
	print "[".localtime(time)."] RUN:  Generating a common alignment statistics file.\n";
	my($header,$output);
	foreach my $newStockName (keys %stocksToAlign) {
		if (-e "$pwd/offspring/$newStockName/$newStockName.alignmentStatistics.tsv") {
			open INF,"$pwd/offspring/$newStockName/$newStockName.alignmentStatistics.tsv";
			my $foo = <INF>;
			$output .= <INF>;
			$header = $foo if !$header;
			close INF;
		} else {
			$output .= "$newStockName\tmissing\n";

		}
	}
	chomp($header);

	open OUTF,">$pwd/out_alignmentStatistics.tsv";
	print OUTF "$header\n";
	print OUTF "$output";
	close OUTF;
}

# Finish up
print "[".localtime(time)."] \n";
print "[".localtime(time)."] DONE: Script completed successfully on the following stocks:\n";
foreach my $newStockName (keys %stocksToAlign) {
	print "[".localtime(time)."] - $newStockName\n";
}
print "[".localtime(time)."] Goodbye.\n";
