#!/usr/bin/perl -w
use strict;
use warnings;    
use File::Spec;
use File::Basename;

use Bio::SearchIO;
use Bio::Seq;
#use Bio::Tools::Run::StandAloneBlast;

# data should be pointed to genome database's root filename generated via getfai.sh
my $data=$ARGV[0];  
my $query=defined($ARGV[1]) ? $ARGV[1] : "./OR_query_mini_ORN_gharial.fasta";
my $query_filename = fileparse($query);
my $output_dir=defined($ARGV[2]) ? $ARGV[2] : '.'; # Default to current directory if not provided
my $db_label = defined($ARGV[3]) ? $ARGV[3] : 'db'; #New change 

# Initialize variables for input and output
#my $out=File::Spec->catfile($output_dir, $query_filename);   $out=~s/\.fasta/\.out/g; # have the output be the same as the query name but replace .fasta with .out
#my $input=$out;
#my $output=$input;
#$output=~s/\.out/\.sum/;

#new change 
my $out = File::Spec->catfile($output_dir, "${db_label}_$query_filename");
$out =~ s/\.(fasta|fa|fna)$/.out/i;

my $input = $out;
my $output = $out;
$input  =~ s/\.out$/.xml/i;   # blast XML that Bio::SearchIO reads
$output =~ s/\.out$/.sum/i;   # your parsed summary table
#new change ends 

my $e="1e-20";
my $pro="tblastn";

# print $query . "\n" . $out . "\n" . $input . "\n" . $output;
# print "$query\n\n$data\n\n$out\n\n";

print "=== BEGIN tblastn ===\n\n";

#new change 2/9
my $cmd = "tblastn -evalue $e -query $query -db $data -out $input -outfmt 5 -num_alignments 200000 -num_descriptions 200000 -num_threads 20";
system($cmd) == 0 or die "BLAST failed: $cmd\n";

#system("tblastn -evalue $e -query $query -db $data -out $input  -outfmt 5 -num_alignments 200000 -num_descriptions 200000 -num_threads 20");
print "=== END tblastn ===\n\n";

#changed $out to $input

  my $file = $input;
  my $in = new Bio::SearchIO(-format =>'blastxml',
                             -file =>$file);
  my $num = $in->result_count;  


  open (OUT, ">>$output");

  print OUT "Query\/Query_length\/Hit\/Hit_length\/E-value\/Bit score\/Percent_identity\/Number_indentity\/Query_Start\/Query_End\/Hit_Start\/Hit_END\/Query_strand\/Hit_strand\n";

     while( my $r = $in->next_result )
   { 
     
     while( my $h = $r->next_hit ) 
     { 
   
         while( my $hsp = $h->next_hsp ) 
                 {
                      my $queryname=$r->query_name;
                       $queryname=~s/\|/\//g;
                       print OUT $queryname,";", " ", $r->query_description,"\/"," ",$r->query_length,"\/";

                       print OUT $h->name, "\/"; 
                       print OUT $hsp->length('total'),"\/", " ",  $hsp->evalue,"\/", $hsp->score, "\/",$hsp->percent_identity ,"\/",$hsp->num_identical,"\/",$hsp->query->start,"\/", " ", $hsp->query->end,"\/"," ", $hsp->hit->start,"\/", " ", $hsp->hit->end,"\/", " ", $hsp->query->strand,"\/", " ",$hsp->hit->strand,"\n";

                 }
                           

      }
}

