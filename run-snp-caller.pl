#!/bin/perl

##########################################################################################
###################################   NOTES  #############################################
##########################################################################################
### This script is designed for finding SNPS in chloroplast genomes and it relies on   ###
### Three applications, picard tools, GATK and VarScan 2 you will need to install these###
### on your system. Then specifiy the path do the jar file for each of these below     ###
### Currently GATK does not work with java 9. If you have a java 9 (or greater) based  ###
### system you will need to specify $java_home below (actually you need to do this even###
### if you have java 8 on your system.                                                 ###
### VarScan2 does not handle duplicated regons very well. As a result you will need to ###
### excise the two inverted repeats from your chloroplast and use this new fine as the ###
### reference_sequence_noirb.fa file. Then use the IRa or ITb recion as the            ###
### IRB_only_ref as input to the script. Also, the script assumes you have paired end  ###
### illumina data and will not work on single end data.                                ###
### Beyond this, as long as you have specified the above, the script will run in an    ###
### automated fashion and will generate the VCF output files for you. However, if you  ###
### want a single VCF output file you will still need to integrate the VCF output for  ###
### the LSC and SSC regions with the VCF output for the two IR regions.                ###
##########################################################################################
##########################################################################################

use Cwd;
use Text::Wrap;
my $cwd = cwd();
$Text::Wrap::columns = 60;

if (!(defined($ARGV[2]))) {
   print "Usage: run-snp-caller.pl <reference_sequence_noirb.fa> <IRB_only_ref> <illumina left reads> <illumina right reads>\n";
   exit;
}

#my $java_home = "java";

#GATK does not work with Java 9, this is just to set an earlier version of java
my $java_home = "/usr/libexec/java_home -v 1.8.0_77 --exec java";


#Set the following 3 to match the location of picard, GATK and VarScan2 on your system
my $picard = "/Users/dyfed/Downloads/picard.jar";
my $gatk = "-Xmx2g -jar /Users/dyfed/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar";
my $varScan = "/Users/dyfed/Downloads/VarScan.v2.3.9.jar";


my $faFil = $ARGV[0];
my $irB = $ARGV[1];
my $lReads = $ARGV[2];
my $rReads = $ARGV[3];
my $sys;

#fix the fasta file just in case.


print "Step 1: Fixing main fasta file\n";
open(INFILE,$faFil) or die $!;

my ($fn,undef) = split /\./, $faFil, 2;
my $nref = $fn."_new.fa";

open(OUTFILE,">${nref}") or die $!;
while (my $line = <INFILE>) {
   chomp($line);
   if ($line =~ /^>/) {
	   $id = $line;
   } else {
	   $seqlin .= $line;
   }
}

print OUTFILE $id . "\n";
print OUTFILE wrap('', '', $seqlin) . "\n";

$faFil = $nref;
close(OUTFILE);
close(INFILE);
$nref = "";
$seqlin = "";

print "Step 1b: Fixing IRB fasta file\n";

open(INFILE,$irB) or die $!;
my ($fn2,undef) = split /\./, $irB, 2;
my $nref = $fn2."_new.fa";

open(OUTFILE,">${nref}") or die $!;
while (my $line = <INFILE>) {
	chomp($line);
	if ($line =~ /^>/) {
		$id = $line;
	} else {
		$seqlin .= $line;
	}
}

print OUTFILE $id . "\n";
print OUTFILE wrap('', '', $seqlin) . "\n";

$irB = $nref;
close(OUTFILE);

#Start with BWA mapping of reads to reference creating a bam file directly with samtools

print "Step 2: Creating BWA index for reference\n";

$sys = "bwa index ".$faFil;
system($sys);

$sys = "bwa index ".$irB;
system($sys);

print "Step 3: Mapping reads to reference\n";

$sys = "bwa mem -M ".$faFil." ".$lReads." ".$rReads." | samtools view -F 4 -Sbh > mapped.bam";
system($sys);
$sys = "bwa mem -M ".$irB." ".$lReads." ".$rReads." | samtools view -F 4 -Sbh > mapped2.bam";
system($sys);


#Add read groups and index (for GATK compatibility)

print "Step 4: Adding read groups for GATK compatibility\n";
$sys = "java -jar ".$picard." AddOrReplaceReadGroups I=mapped.bam O=M7_mapped_rg.bam RGID=ga1 RGLB=LCP RGPL=illumina RGPU=id1 RGSM=1";
system($sys);

$sys = "java -jar ".$picard." AddOrReplaceReadGroups I=mapped2.bam O=M7_mapped2_rg.bam RGID=ga1 RGLB=LCP RGPL=illumina RGPU=id1 RGSM=1";
system($sys);


print "Step 5: sorting and indexing the mapping BAM file\n";
$sys = "samtools sort M7_mapped_rg.bam > M7_sorted.bam";
system($sys);

$sys = "samtools sort M7_mapped2_rg.bam > M7_sorted2.bam";
system($sys);

$sys = "samtools index M7_sorted.bam";
system($sys);

$sys = "samtools index M7_sorted2.bam";
system($sys);



#mark and remove duplicates (picard tools)

print "Step 6: marking and removing optical duplicates\n";
$sys = "java -jar ".$picard." MarkDuplicates I=M7_sorted.bam O=M7_deduped.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT";
system($sys);

$sys = "java -jar ".$picard." MarkDuplicates I=M7_sorted2.bam O=M7_deduped2.bam M=marked_dup_metrics2.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT";
system($sys);

$sys = "samtools sort M7_deduped.bam > M7_dedup_sorted.bam";
system($sys);

$sys = "samtools sort M7_deduped2.bam > M7_dedup2_sorted.bam";
system($sys);

$sys = "samtools index M7_dedup_sorted.bam";
system($sys);

$sys = "samtools index M7_dedup2_sorted.bam";
system($sys);


#re-align indels
print "Step 7: re-aligning indel regions with GATK\n";
$sys = "samtools faidx ".$faFil;
system($sys);

$sys = "samtools faidx ".$irB;
system($sys);

$sys = "java -jar ".$picard." CreateSequenceDictionary R=".$faFil." O=".$fn."_new.dict";
system($sys);

$sys = "java -jar ".$picard." CreateSequenceDictionary R=".$irB." O=".$fn2."_new.dict";
system($sys);

$sys = $java_home." ".$gatk." -T RealignerTargetCreator -R ".$faFil." -I M7_dedup_sorted.bam -o M7_out.intervals --filter_reads_with_N_cigar";
system($sys);

$sys = $java_home." ".$gatk." -T RealignerTargetCreator -R ".$irB." -I M7_dedup2_sorted.bam -o M7_out2.intervals --filter_reads_with_N_cigar";
system($sys);

$sys = $java_home." ".$gatk." -I M7_dedup_sorted.bam -R ".$faFil." -T IndelRealigner -targetIntervals M7_out.intervals -o M7_realigned.bam --filter_reads_with_N_cigar";
system($sys);

$sys = $java_home." ".$gatk." -I M7_dedup2_sorted.bam -R ".$irB." -T IndelRealigner -targetIntervals M7_out2.intervals -o M7_realigned2.bam --filter_reads_with_N_cigar";
system($sys);

print "\n";
print "Step 8: calling SNPs with VarScan2\n";
print "Calling SNPs on main sequence\n";
$sys = "samtools mpileup -B -f ".$faFil." M7_realigned.bam | java -jar ".$varScan." mpileup2cns --min-coverage 8 --min-avg-qual 20 --min-var-freq 0.2 --p-value 0.005 --variants --strand-filter 0 --min-strand2 0 --output-vcf > variants.vcf";
system($sys);

print "\n";
print "Calling SNPs on IRb sequence\n";

$sys = "samtools mpileup -B -f ".$irB." M7_realigned2.bam | java -jar ".$varScan." mpileup2cns --min-coverage 8 --min-avg-qual 20 --min-var-freq 0.2 --p-value 0.005 --variants --strand-filter 0 --min-strand2 0 --output-vcf > variants2.vcf";
system($sys);

print "Pipeline completed\n";

#samtools mpileup -B -f Saccharum_officinarum_BH10_12.fa BH10_12_realigned.bam | java -jar /Users/dyfed/Downloads/VarScan.v2.3.9.jar mpileup2cns --min-coverage 10 --min-avg-qual 20 --min-var-freq 0.5 --p-value 0.005 --variants --output-vcf > variants.vcf


