#!/bin/perl

use Cwd;
use Text::Wrap;
use Bio::SeqIO;

my $cwd = cwd();
$Text::Wrap::columns = 60;

my %idHash;
#my $bwa = "/home/dyfed/Downloads/bwa.kit/bwa";
#my $samtools = "/home/dyfed/Downloads/samtools-1.3.1/samtools";
#my $fqdump = "/home/dyfed/Downloads/sratoolkit.2.8.0-centos_linux64/bin/fastq-dump";
my $bwa = "bwa";
my $samtools = "samtools";
my $fqdump = "fastq-dump";
my $ctr = 0;
my $maxnum;
my $spec = undef;
my $comm;

if (!(defined($ARGV[0]))) {
   print "Usage: build-genes.pl <ref> <sras>\n";
   exit;
}

my $rroot = $ARGV[0];

open(INFILE,$rroot) or die $!;

   while (my $line = <INFILE>) {
   chomp($line);
   $ctr++;
   if (!($line =~ /^\t/)) {
	   if ($ctr == 1) {
#	   if ((!(defined($spec))) && (!($spec eq ""))) {
	      $spec = $line;
	      print "Species $spec \n";
	      $spec =~ s/ /_/g;
           } else {
		   &process_sra($spec);
		   $spec = $line;
		   $spec =~ s/ /_/g;
	   }
   } else {
	   $sra = $line;
	   $sra =~ s/\t//;
	   my $sstr = substr $sra, 0, 6;
	   if ($sra =~ /^SRR/) {
	      $comm = "ncftpget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR";
           } elsif ($sra =~ /^ERR/) {
	      $comm = "ncftpget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR";
           } elsif ($sra =~ /^DRR/) {
	      $comm = "ncftpget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/DRR";
           }
      $comm .= "/".$sstr."/".$sra."/".$sra.".sra";

      system($comm);
      system($comm);
      system($comm);
      print $comm . "\n";

      }

}     

								 






#unpack the SRA






sub process_sra{
    my $species = shift @_;

#Fix the fasta file and header

      my $ref = $species.".fa";

      my $seqio = Bio::SeqIO->new(-file => $ref, -format => "fasta");
      while(my $seq = $seqio->next_seq) {
		       #throw error if not fasta file
      }

      my $nref = $species."_new.fa";

      open(INFILE,$ref) or die $!;
      open(OUTFILE,">${nref}") or die $!;
      print "Fixed reference file: $nref \n";

      while (my $line = <INFILE>) {
	  chomp($line);
	  if ($line =~ /^>/) {
		$id = $line;
	} else {
	        $seqlin .= $line;
	}

     }

    print OUTFILE ">" .$species . "\n";
    print OUTFILE wrap('', '', $seqlin) . "\n";
    close(OUTFILE);

    my $sys = $bwa." index ".$nref;
    system($sys);

    opendir(DH, $cwd);
    my @files = readdir(DH);
    closedir(DH);

    foreach my $file (@files) {
       chomp($file);
       if ($file =~ /\.sra$/) {
          my ($s1,undef) = split /\./, $file, 2;
          $idHash{$s1} = $s1;
       }
    }


my $ctr = 0;

my $bamlist;

foreach my $keys (keys(%idHash)) {
  $ctr++;
  $sys = $fqdump." --split-files ".$keys.".sra";
  print $sys . "\n";

  system($sys);

  my $file = $keys."_2.fastq";

  if (-f $file) {
     $sys = $bwa." mem -M ".$nref." ".$keys."_1.fastq ".$keys."_2.fastq > mapped".$ctr.".sam";
  } else {
     $sys = $bwa." mem -M ".$nref." ".$keys."_1.fastq > mapped".$ctr.".sam";
  }

  print "BWAMAP $sys \n";
  system($sys);

  $sys = $samtools." view -F 4 -Sbh mapped".$ctr.".sam > mapped".$ctr.".bam";
  system($sys);

  $sys = "rm mapped".$ctr.".sam";
  system($sys);

  $bamlist .= " mapped".$ctr.".bam";
}

$sys = $samtools." merge merged.bam ".$bamlist;
print "Merging: ".$sys."\n";
system($sys);

$sys = $samtools." sort merged.bam > merged_sorted.bam";
system($sys);

$sys = $samtools." index merged_sorted.bam";
system($sys);

$sys = $samtools." depth -a merged_sorted.bam > merged_sorted.count";
system($sys);

#prep the count file for R
open(INFILE,"<merged_sorted.count");
my $nref = "merged_sorted_new.cov";

open(OUTFILE,">${nref}") or die $!;

print OUTFILE "CHR LOC LOG10\n";

while (my $line = <INFILE>) {
   chomp($line);
   $ctr++;
   my ($chr,$loc,$cov) = split /\t/, $line, 3;
   if (!($cov == 0)) {
      $cov = log($cov)/log(10);
   } else {$cov = 0;}

   print OUTFILE $chr." ".$loc." ".$cov."\n";
}

#Perform cleanup
$sys = "mkdir ".$species;
system($sys);

$sys = "rm *.sra";
system($sys);

$sys = "gzip *.fastq";
system($sys);

$sys = "mv *.bam *.idx *.count *.cov ".$spcies."*.fa.* ".$species;
system($sys);

$sys = "mv *.fastq.gz ".$species;
system($sys);

}

1;
