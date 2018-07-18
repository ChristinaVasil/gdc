#!/bin/env perl
 
use strict;
use warnings;


open INPUT, "rnaseqNormalized.txt" or die "cannot open file";
open OUTPUT, ">rnaseqNormlog.txt" or die "cannot open file";


my @row;
my @newrow;
my $firstrow;

while($b=<INPUT>){

chomp $b;

if ($. == 1){

	$firstrow = $b;
	print OUTPUT "$firstrow\n";
	
}

if ($. > 1){
	

	@row= split ("\t", $b);
	 my $length= @row;
	for(my $i=1; $i<$length; $i++){
		
		$newrow[$i-1] = $row[$i];
		
	}
	
	
	my @logged= log2(@newrow);
	my $mergedlogrow= join("\t", @logged);
	print OUTPUT "$row[0]\t$mergedlogrow\n";
}
}


close INPUT;
close OUTPUT;

exit(0);
 

sub log2 {
   my @n = @_; 
   my @t = ();
   foreach my $n (@n){
      if ($n == 0){
         $n = '0.5';
      }
      my $t = log($n)/log(2);
      #rounded to two decimal places
     # $t = sprintf("%.2f",$t);
      push(@t,$t);
   }
   return(@t);
}