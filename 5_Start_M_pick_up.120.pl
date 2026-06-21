#!/usr/bin/perl -w

use File::Basename;
use File::Spec;

# Check if there are any arguments provided
if (@ARGV == 0) {
    die "***************************** ERROR *****************************\nERROR: Please provide absolute path of spp directory\n";
}

my $SPP_PATH = $ARGV[0];
my $first_tm = $ARGV[1]; 

# Check if the file exists
unless (-d $SPP_PATH) {
    die "***************************** ERROR *****************************\nERROR: File does not exist at '$SPP_PATH'.\n";
}

# Check if the second argument is an integer
if ($first_tm !~ /^\d+$/) {
    # If it's not an integer, try to convert it to an integer
    if ($first_tm =~ /^(\d+)$/) {
        $first_tm = $1;
    } else {
        die "The second argument '$first_tm' is not a valid integer.\n";
    }
}

my $spp = basename($SPP_PATH);

print "'$SPP_PATH'\n'$spp'";

opendir(FOLDER,"$SPP_PATH") or die print "can not find the file";
my @array = grep(/${spp}_step_4_result_mafftAlignment_removed.fasta$/,readdir(FOLDER)); 
close FOLDER;


foreach my $filename ( @array ){

    my $input=File::Spec->catfile($SPP_PATH, $filename);
    my $fasta_output= File::Spec->catfile($SPP_PATH, "step_5_result.pickedM.fasta");
    my $text_output= File::Spec->catfile($SPP_PATH, "step_5_result.pickedM.CoordinatesOfBestStart.txt");

    open (FILE,"$input");
    open (FASTA, ">$fasta_output");
    open (TEXT, ">$text_output");

    my @array; 
    my $species;
    my %sequences;

    while( <FILE> )	{
        my $line=$_;  chomp $line;
        if( $line =~ /^>(.+)/ ){
            $species = $1;
            $sequences{$species} = '';  push(@array, $species) ;
        }
        else{
            $sequences{$species} .= $line;                      
        }

    }
    close FILE;
	
	foreach my $name( @array){ 
	    my $lengthall=length($sequences{$name});
         
		my $part2=substr($sequences{$name}, ($first_tm-1), ($lengthall-$first_tm+1));
           $part2=~s/-//g;		
	
	    my $Nsite=substr($sequences{$name}, 0, ($first_tm-1));
            $Nsite=~s/-//g;
            $Nsite=~s/\n//g;
		
		if ($Nsite){
		
            my $length=length($Nsite);
            my @array35; my @array2034; my @array21;
            my $loc=0;
                    
            for (my $i=($length-1); $i>=0; $i--){
                $loc++;
                my $aa=substr($Nsite, $i, 1); # print OUT $aa,"\n";
                if ($aa=~/M/){
                    if (($loc>20) and ($loc<35)) {push (@array2034,$loc)}
                    if ($loc>=35) {push (@array35, $loc)}
                    if ($loc<=20) {push (@array21, $loc)}
                }
            }
		
            my $best_start;

            print FASTA ">$name\n";
            print TEXT ">$name\t";	

            if (@array2034)  {
                $best_start=min(@array2034);
                my $Nsequence=substr($Nsite, ($length-$best_start), $best_start);

                print FASTA $Nsequence;
                print FASTA $part2,"\n\n";

                print TEXT ($length-$best_start), "\n";
                
                next;
            }
            
            elsif (@array35)  {
                $best_start=min(@array35);
                my $Nsequence=substr($Nsite, ($length-$best_start), $best_start);
                
                print FASTA $Nsequence;
                print FASTA $part2,"\n\n";

                print TEXT ($length-$best_start), "\n";
                
                next;
            }
		
            elsif (@array21)  {
                $best_start=max(@array21);
                my $Nsequence=substr($Nsite, ($length-$best_start), $best_start);
                
                print FASTA $Nsequence;
                print FASTA $part2,"\n\n";

                print TEXT ($length-$best_start), "\n";

                next;
            }
		}
	}
}
	
	
	
		  
	  sub max {
   my($max_so_far) = shift @_;
   foreach (@_){                         
      if($_>$max_so_far){                  
           $max_so_far=$_;
      }
   }
   return $max_so_far;                      
}

   

sub min {
   my($min_so_far) = shift @_;
   foreach (@_){                         
      if($_<$min_so_far){                  
           $min_so_far=$_;
      }
   }
   return $min_so_far;                      
}
		
