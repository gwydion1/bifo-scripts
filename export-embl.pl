#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use IO::String;
use Bio::SeqIO;
use DateTime;
use Capture::Tiny ':all';


use Text::Wrap qw(wrap);
$Text::Wrap::unexpand = 0;
my $cwd = cwd();

if (!(defined($ARGV[1]))) {
   print "Usage: export-embl.pl <tbl-file> <extras-file> <classification> <fastafil>\n";
   exit;
}

my $tblfil = $ARGV[0];
my $extras = $ARGV[1];
my $class = $ARGV[2];
my $fastafil = $ARGV[3];
my $dna = "";
my $id;
my $found = undef;
my $ctr = 0;
my %geneHash;
my %sortHash;
my $gname;
my $val;
my %nameHash;
my $seq;
my $slen;
my $embseq;
my $gtctr = 0;
my $gline;
my $tline;
my $cdsline;
my $rplctr = 0;

open(INFILE, "gene-names.in");

while (my $line = <INFILE>)
{
   chomp($line);
   my ($gid,$dscs) = split /\t/, $line, 2;
   $nameHash{$gid} = $dscs;
}

close(INFILE);

open(INFILE, "<${fastafil}");

while (my $line = <INFILE>) {
   chomp($line);
   if (!($line =~ /^>/)) {
      $seq .= $line;
      $slen = $slen + length($line);
   }
}


open(INFILE, "<${tblfil}");
my ($fnm,undef) = split /\./, $tblfil, 2;
my $embl = $fnm."\.embl";

open(OUTFILE, ">${embl}");


my $of2 = $fnm."\.chr";
open(OUTFILE2, ">${of2}");

print OUTFILE2 $fnm."-pt I chromosome pastid\n";


      while (my $line = <INFILE>)
      {
      chomp($line);
      if (!($line =~ /~/)) {
	 if (!($line =~ /intron/)) {
            print $line . "\n";
            my ($iid,$gst,$gen,$ori) = split /\t/, $line, 4;
            my ($id,undef) = split /_exon/, $iid, 2;
	    if ($id =~ /3end/) {
               my ($nid,undef) = split /_3end/, $id, 2;
	       $id = $nid;
            }
	    if ($ori eq '+') {
	       if ($id =~ /rpl23/) {
	          if ($rplctr == 0) {
		     $id .= "1";
		     $rplctr = 10;
		  }
	       }
	       $gname = $id . "-plus";
            } else {
               $gname = $id . "-minus";
	    }
	    if ($id eq "rps12_5end") {
	       $gname = "rps12-plus";
	       $geneHash{$gname} = $gst."^".$gen."^".$ori;
	       $gname = "rps12-minus";
	    }
	    print $gname . "\t" . $line . "\n";
	    if (!(exists $geneHash{$gname})) {
               $geneHash{$gname} = $gst."^".$gen."^".$ori;
            } else {
	       my $ent = $geneHash{$gname};
	       $ent .= "::". $gst."^".$gen."^".$ori;
	       $geneHash{$gname} = $ent;
	       $ent = undef;
	   }
      }  
      }

}

close(INFILE);

open(INFILE, "<${extras}");


while (my $line = <INFILE>) {
   chomp($line);
   my ($iid,$gst,$gen,$ori) = split /\t/, $line, 4;
   if ($iid =~ /^ORF/) {
      $iid = lc $iid;
   }
   if ($iid =~ /^orf68/) {
	   $iid = "orf69";
   }
   if ($iid =~ /^ycf/) {
	   $iid =~ s/_/ /;
   }
   if ($ori eq '+') {
       $gname = $iid . "-plus";
   } else {
       $gname = $iid . "-minus";
   }
   if (!($gname =~ /rpl23/)) {
      if (!(exists $geneHash{$gname})) {
         $geneHash{$gname} = $gst."^".$gen."^".$ori;
      } else {
         my $ent = $geneHash{$gname};
         $ent .= "::". $gst."^".$gen."^".$ori;
         $geneHash{$gname} = $ent;
         $ent = undef;
      }
   }
}


close(INFILE);

my $species;
my $cv;
my $sv;
my $oc;
my $tax;
my $project;


open(INFILE, "<${class}");

while (my $line = <INFILE>) {
   chomp($line);
   if ($line =~ /^species/) {
      (undef,$species) = split /:/, $line, 2;
   } elsif ($line =~ /^cv:/) {
      (undef,$cv) = split /:/, $line, 2;
   } elsif ($line =~ /^sv:/) {
      (undef,$sv) = split /:/, $line, 2;
   }  elsif ($line =~ /^oc:/) { 
      (undef,$oc) = split /:/, $line, 2;
      $oc =~ s/^://;
   } elsif ($line =~ /^tax:/) {
      (undef,$tax) = split /:/, $line, 2;
   } elsif ($line =~ /^project:/) {
      (undef,$project) = split /:/, $line, 2;
   }
}

foreach my $key (keys %geneHash) {
#	print $key . "\t" . $geneHash{$key} . "\n";
	my ($start,undef) = split /\^/, $geneHash{$key}, 2;
	$val = $key."^^".$geneHash{$key};
	if (!(exists($sortHash{$start}))) {
           $sortHash{$start} = $val;
        } else {
	   $val = $sortHash{$start};
           $val .= "**" . $key."^^".$geneHash{$key};
	   $sortHash{$start} = $val;
        
       }
}


#start printing the header stuff

my $dt = DateTime->now;
my $pdt = $dt->day."-".$dt->month_abbr."-".$dt->year;

print OUTFILE "ID   XXX; XXX; circular; XXX; XXX; XXX; XXX.\n";
print OUTFILE "XX\n";
print OUTFILE "AC   ".$fnm."-pt;\n";
print OUTFILE "XX\n";
if ((defined($project)) && (!($project eq ''))) {
	print OUTFILE "PR   Project:".$project.";\n";
} else {
   print OUTFILE "PR   Project:PRJEB1234;\n";
}
print OUTFILE "XX\n";
if ((defined($cv)) && (!($cv eq ''))) {
   print OUTFILE "DE   ".$species." chloroplast complete genome, cultivar ".$cv."\n";
} elsif ((defined($sv)) && (!($sv eq ''))) {
   print OUTFILE "DE   ".$species." chloroplast complete genome, specimen voucher ".$sv."\n";	
} else {
   print OUTFILE "DE   ".$species." chloroplast complete genome\n";
}
print OUTFILE "XX\n";
print OUTFILE "KW   complete genome.\n";
print OUTFILE "XX\n";
print OUTFILE "OS   ".$species."\n";
my $max = 71;
local $Text::Wrap::columns = $max;
print OUTFILE wrap('OC   ', 'OC   ', $oc);
print OUTFILE "\n";
print OUTFILE "XX\n";
print OUTFILE "RN   [1]\n";
print OUTFILE "RP   1-".$slen."\n";
print OUTFILE "RA   Lloyd Evans D.;\n";
print OUTFILE "RT   ;\n";
print OUTFILE "RL   Submitted (".$pdt.") to INSDC.\n";
print OUTFILE "RL   South African Sugarcane Research Institute, 170 Flanders Drive, Mount\n";
print OUTFILE "RL   Edgecombe, Durban 4300, South Africa,\n";
print OUTFILE "XX\n";
print OUTFILE "CC   ##Assembly-Data-START##\n";
print OUTFILE "CC   Assembly Method       :: SPADES v. 3.9\n";
print OUTFILE "CC   Finishing Method      :: Pilon\n";
print OUTFILE "CC   Sequencing Technology :: Illumina\n";
print OUTFILE "CC   ##Assembly-Data-END##\n";
print OUTFILE "XX\n";
print OUTFILE "FH   Key             Location/Qualifiers\n";
print OUTFILE "FH\n";
print OUTFILE "FT   source          1..".$slen."\n";
if ((defined($cv)) && (!($cv eq ''))) {
      print OUTFILE "FT                   /organism=\"".$species." cultivar ".$cv."\"\n";
} else {
      print OUTFILE "FT                   /organism=\"".$species."\"\n";
}
print OUTFILE "FT                   /organelle=\"plastid:chloroplast\"\n";
if ((defined($cv)) && (!($cv eq ''))) {
	print OUTFILE "FT                   /cultivar=\"".$cv."\"\n";
}
if ((defined($sv)) && (!($sv eq ''))) {
	print OUTFILE "FT                   /specimen_voucher=\"".$sv."\"\n";
}
print OUTFILE "FT                   /mol_type=\"genomic DNA\"\n";
if ((defined($tax)) && (!($tax eq ''))) {
	print OUTFILE "FT                   /db_xref=\"taxon:".$tax."\"\n";
}
print OUTFILE "FT   misc_feature    complement(1..".$slen.")\n";
print OUTFILE "FT                   /note=\"JLA, junction IRA-LSC\"\n";



for my $k (sort {$a<=>$b} keys %sortHash) {
	  print $k . "\t" . $sortHash{$k}. "\n";
	  my @lines = split /\*\*/, $sortHash{$k};

	  foreach my $line (@lines) {
	     my ($gid,$exes) = split /\^\^/, $line, 2;
	     my @exons = split /::/, $exes;
	     my ($gnnm,undef) = split /-minus/, $gid, 2;
	     my ($gnm,undef) = split /-plus/, $gnnm, 2;
	     my $gdsc = $nameHash{$gnm};
	     print "Gene: " . $gnm . "\t" . $exes . "\t" . $gdsc . "\n";
	     if ((scalar @exons) == 1) {
		     my $exon = shift(@exons);
                my ($fst,$fen,$fori) = split /\^/, $exon;
		if ($gnm eq "rpl231") {
                   $gnm = "rpl23";
		   $gdsc = $nameHash{$gnm};
	        }
	        if ($gnm eq "LSC") {
                   print OUTFILE "FT   misc_feature    ".$fst."..".$fen."\n";
		   print OUTFILE 'FT                   /note="large single copy region; LSC"'."\n";
	        } elsif ($gnm eq "FULL") {
                   #do nothing
	        } elsif ($gnm eq "IRB") {
		   my $nlc = $fst -1;
		   print OUTFILE "FT   misc_feature    ".$nlc."..".$fst."\n";
		   print OUTFILE 'FT                   /note="JLB, junction LSC-IRB"'."\n";
		   print OUTFILE "FT   misc_feature    ".$fst."..".$fen."\n";
		   print OUTFILE 'FT                   /note="IRB, Inverted Repeat Region, Copy B"'."\n"
	       } elsif ($gnm eq "IRA") {
		   my $nlc = $fst -1;
		   print OUTFILE "FT   misc_feature    ".$nlc."..".$fst."\n";
		   print OUTFILE 'FT                   /note="JSA, junction SSC-IRA"'."\n";
	           print OUTFILE "FT   misc_feature    ".$fst."..".$fen."\n";
		   print OUTFILE 'FT                   /note="IRA, Inverted Repeat Region, Copy A"'."\n"
	       } elsif ($gnm eq "SSC") {
		   my $nlc = $fst -1;
		   print OUTFILE "FT   misc_feature    ".$nlc."..".$fst."\n";
		   print OUTFILE 'FT                   /note="JSB, junction IRB-SSC"'."\n";
		   print OUTFILE "FT   misc_feature    ".$fst."..".$fen."\n";
		   print OUTFILE 'FT                   /note="SSC, short single copy region; SSC"'."\n"
	      } elsif ($gnm =~ /^oriA/) {
                   print OUTFILE "FT   rep_origin      ".$fst."..".$fen."\n";
		   print OUTFILE "FT                   /note=\"".$gnm."\"\n";
	      } elsif ($gnm =~ /pseudo/) {
	         if ($gnm =~ /^ACRS/) {
		    if ($gnnm =~ /-plus/) {
                       print OUTFILE "FT   gene            complement(<".$fst."..".$fen.")\n";
		       print OUTFILE 'FT                   /gene="ACRS"'."\n";
		       print OUTFILE "FT                   /pseudo\n";
                       print OUTFILE "FT   CDS             complement(<".$fst."..".$fen.")\n";
		       print OUTFILE 'FT                   /gene="ACRS"'."\n";
		       print OUTFILE "FT                   /codon_start=1\n";
                       print OUTFILE "FT                   /pseudo\n";
		       print OUTFILE "FT                   /note=\"ACR toxin sensitivity protein pseudogene gained from the mitochondrion\"\n";
	            } else {
                       print OUTFILE "FT   gene            <".$fst."..".$fen."\n";
		       print OUTFILE 'FT                   /gene="ACRS"'."\n";
		       print OUTFILE "FT                   /pseudo\n";
		       print OUTFILE "FT   CDS             <".$fst."..".$fen."\n";
		       print OUTFILE 'FT                   /gene="ACRS"'."\n";
		       print OUTFILE "FT                   /codon_start=1\n";
		       print OUTFILE "FT                   /pseudo\n";
		       print OUTFILE "FT                   /note=\"ACR toxin sensitivity protein pseudogene gained from the mitochondrion\"\n";
		    }
	         } elsif ($gnm =~ /^rpl2-p/) {
		    print OUTFILE "FT   gene            complement(<".$fst."..>".$fen.")\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT   CDS             complement(<".$fst."..>".$fen.")\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE "FT                   /codon_start=2\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT                   /note=\"pseudogene for ribosomal protein L2\"\n";
		 } elsif ($gnm =~ /^tRNA-Thr/) {
                    print OUTFILE "FT   gene            <".$fst."..>".$fen."\n";
		    print OUTFILE 'FT                   /gene="tRNA-Thr"'."\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT   tRNA            <".$fst."..>".$fen."\n";
		    print OUTFILE 'FT                   /gene="tRNA-Thr"'."\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT                   /note=\"trnT pseudogene for tRNA-Thr(GGU)\"\n";
		 } elsif ($gnm =~ /^rpl23-p/) {
                    print OUTFILE "FT   gene            <".$fst."..".$fen."\n";
		    print OUTFILE 'FT                   /gene="rpl23"'."\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT   CDS             <".$fst."..".$fen."\n";
		    print OUTFILE 'FT                   /gene="rpl23"'."\n";
		    print OUTFILE "FT                   /codon_start=1\n";
		    print OUTFILE "FT                   /pseudo\n";
		    print OUTFILE "FT                   /note=\"rpl23 pseudogene for ribosomal protein L23\"\n";
		 }
	       } else {
                if ($gnm =~ /^trn/) {
                   if ($fori eq "-") {
                      print OUTFILE "FT   gene            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   tRNA            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$gdsc.'"'."\n";
		   } else {
                      print OUTFILE "FT   gene            ".$fst."..".$fen."\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   tRNA            ".$fst."..".$fen."\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$gdsc.'"'."\n";
		   }
		} elsif ($gnm =~ /^rrn/) {
		      if ($gnm eq "rrn45") {
                         $gnm = 'rrn.4.5';
			 $gdsc = '4.5s ribosomal RNA';
		      }
                   if ($fori eq "-") {
		      print OUTFILE "FT   gene            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   rRNA            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$gdsc.'"'."\n";
		   } else {
		      print OUTFILE "FT   gene            ".$fst."..".$fen."\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   rRNA            ".$fst."..".$fen."\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$gdsc.'"'."\n";
		  }



		} else {
# add gene code here
                   my ($prod,$note) = split /\t/, $gdsc, 2;
		   if ($fori eq "-") {
	              print OUTFILE "FT   gene            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   CDS             complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$prod.'"'."\n";
		      if ((defined $note) && (!($note eq ''))) {
			 $note = '/note="'.$note.'"'."\n";
			 my $max = 80;
			 local $Text::Wrap::columns = $max;
			 print OUTFILE wrap('FT                   ', 'FT                   ', $note);
			 # print OUTFILE 'FT                   /note="'.$note.'"'."\n";
		      }
		      print OUTFILE "FT                   /codon_start=1\n";
		      print OUTFILE "FT                   /transl_table=11\n";
		   } else {
                      print OUTFILE "FT   gene            complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE "FT   CDS             complement(".$fst."..".$fen.")\n";
		      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		      print OUTFILE 'FT                   /product="'.$prod.'"'."\n";
# Apply Text::Wrap or something to fix long note descriptors.
		      if ((defined $note) && (!($note eq ''))) {
			 $note = '/note="'.$note.'"'."\n";
			 my $max = 80;
			 local $Text::Wrap::columns = $max;
			 print OUTFILE wrap('FT                   ', 'FT                   ', $note);

#			 print OUTFILE 'FT                   /note="'.$note.'"'."\n";
		      }
		      print OUTFILE "FT                   /codon_start=1\n";
		      print OUTFILE "FT                   /transl_table=11\n";

		   }
		}


	      }
	   } else {

              #process multi-exon genes
              if ($gid =~ "rps12") {
                 my $eex1 = shift(@exons);
		 my $eex2 = shift(@exons);
		 my $eex3 = shift(@exons);
		 my ($e1st,$e1en,undef) = split /\^/, $eex1;
		 my ($e2st,$e2en,undef) = split /\^/, $eex2;
		 my ($e3st,$e3en,undef) = split /\^/, $eex3;
		 if ($gid =~ /-minus/) {
		   print OUTFILE "FT   gene            complement(join(".$e2st."..".$e3en.",".$e1st."..".$e1en."))\n";
		   print OUTFILE 'FT                   /gene="rps12"'."\n";
		   print OUTFILE "FT   CDS             complement(join(".$e2st."..".$e2en.",".$e3st."..".$e3en.",".$e1st."..".$e1en."))\n";
		   print OUTFILE "FT                   /codon_start=1\n";
		   print OUTFILE "FT                   /transl_table=11\n";
		   print OUTFILE "FT                   /trans_splicing\n";
		   print OUTFILE 'FT                   /gene="rps12"'."\n";
		   print OUTFILE 'FT                   /product="ribosomal protein S12"'."\n";
		 } elsif ($gid =~ /-plus/) {
                   print OUTFILE "FT   gene            join(complement(".$e1st."..".$e1en."),".$e2st."..".$e3en.")\n";
                   print OUTFILE 'FT                   /gene="rps12"'."\n";
		   print OUTFILE "FT   CDS             join(complement(".$e1st."..".$e1en."),".$e2st."..".$e2en.",".$e3st."..".$e3en.")\n";
		   print OUTFILE "FT                   /codon_start=1\n";
		   print OUTFILE "FT                   /transl_table=11\n";
		   print OUTFILE "FT                   /trans_splicing\n";
		   print OUTFILE 'FT                   /gene="rps12"'."\n";
		   print OUTFILE 'FT                   /product="ribosomal protein S12"'."\n";
		 }


	      } elsif ($gnm eq "rpl2") {
                 my $eex1 = shift(@exons);
		 my $eex2 = shift(@exons);
		 my ($e1st,$e1en,undef) = split /\^/, $eex1;
		 my ($e2st,$e2en,undef) = split /\^/, $eex2;
		 if ($gid =~ /-minus/) {
                    my $mst = $e2en -2;
		    print OUTFILE "FT   gene            complement(".$e1st."..".$e2en.")\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE "FT   CDS             complement(join(".$e1st."..".$e1en.",".$e2st."..".$e2en."))\n";
		    print OUTFILE "FT                   /codon_start=1\n";
		    print OUTFILE "FT                   /transl_table=11\n";
		    print OUTFILE "FT                   /transl_except=(pos:".$mst."..".$e2en.",aa:Met)\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE 'FT                   /product="ribosomal protein L2"'."\n";
		    print OUTFILE 'FT                   /note="50S ribosomal protein L2"'."\n";
		    print OUTFILE 'FT                   /note="acg start"'."\n";
		 } elsif ($gid =~ /-plus/) {
                    my $men = $e1st +2;
		    print OUTFILE "FT   gene            ".$e1st."..".$e2en."\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE "FT   CDS             join(".$e1st."..".$e1en.",".$e2st."..".$e2en.")\n";
		    print OUTFILE "FT                   /codon_start=1\n";
		    print OUTFILE "FT                   /transl_table=11\n";
		    print OUTFILE "FT                   /transl_except=(pos:".$e1st."..".$men.",aa:Met)\n";
		    print OUTFILE 'FT                   /gene="rpl2"'."\n";
		    print OUTFILE 'FT                   /product="ribosomal protein L2"'."\n";
		    print OUTFILE 'FT                   /note="50S ribosomal protein L2"'."\n";
		    print OUTFILE 'FT                   /note="acg start"'."\n";
		 }

	      } else {
		      #process all other genes
                      my $ctr = 0;
		      my $gst;
		      my $gen;
		      my $exstr;
		      my ($prod,$note) = split /\t/, $gdsc, 2;
		      foreach my $exon (@exons) {
			      $ctr++;
			      my ($exst,$exen,$ori) = split /\^/, $exon, 3;
			      $gen = $exen;
			      $exstr .= $exst . ".." . $exen . ",";
			      if ($ctr == 1) {
				      $gst = $exst;
			      }
		      }
		      $exstr =~ s/,+$//;
		      if ($gid =~ /-minus/) {
			      print OUTFILE "FT   gene            complement(".$gst."..".$gen.")\n";
			      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
			      if ($gid =~ /^trn/) {
				 print OUTFILE "FT   tRNA            complement(join(".$exstr."))\n";
			      } else {
			         print OUTFILE "FT   CDS             complement(join(".$exstr."))\n";
		              
			         print OUTFILE "FT                   /codon_start=1\n";
			         print OUTFILE "FT                   /transl_table=11\n";
			      }
			      print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
			      print OUTFILE 'FT                   /product="'.$prod.'"'."\n";

		     } elsif ($gid =~ /-plus/) {
                        print OUTFILE "FT   gene            ".$gst."..".$gen."\n";
		        print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
			if ($gid =~ /^trn/) {
			   print OUTFILE "FT   tRNA            join(".$exstr.")\n";
			} else {
		           print OUTFILE "FT   CDS             join(".$exstr.")\n";
		           print OUTFILE "FT                   /codon_start=1\n";
		           print OUTFILE "FT                   /transl_table=11\n";
		       }
		        print OUTFILE 'FT                   /gene="'.$gnm.'"'."\n";
		        print OUTFILE 'FT                   /product="'.$prod.'"'."\n";	

		     }

                     if ((defined $note) && (!($note eq ''))) {
	                $note = '/note="'.$note.'"'."\n";
		        my $max = 80;
			local $Text::Wrap::columns = $max;
			print OUTFILE wrap('FT                   ', 'FT                   ', $note);

#                        print OUTFILE 'FT                   /note="'.$note.'"'."\n";
                    }

	     }	      

	  }


	     $gdsc = undef;

	     
          }
  }


  #process sequence
#print $seq . "\n";


my $in = Bio::SeqIO->new(-string => $seq);

my $ino = $in->next_seq();

my $eseq;
my $stringfh = IO::String->new($eseq);
my $outObj = new Bio::SeqIO(-format => 'EMBL', -fh  => $stringfh);
$outObj->write_seq($ino);

open my $fh, '<', \$eseq or die $!;
while (my $seqline = <$fh>) {
   if ($seqline =~ /^ / || $seqline =~ /^SQ / || $seqline =~ /\/\//) {
      $embseq .= $seqline;
   }
}

print OUTFILE "XX\n";
print OUTFILE $embseq;


1;

