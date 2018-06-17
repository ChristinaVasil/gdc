#!usr/bin/perl -w

use strict;
use warnings;
use Set::Object;
use Set::Object qw(set);

# open mapping file
open INPUT, "gdc_sample_sheet.2018-06-06.tsv" or die "Cannot find file";

# read case ids

my @case_id=();
my $row_counter=0; 
# initialize an empty hash patient_file_association, value is case id and key is the associated file which is unique 
my %patient_methylation_file_association;

while (my $line = <INPUT>){
	
	if ($. > 1){
	
		my @lines = split ("\t", $line);
		
		if($lines[2] eq 'DNA Methylation'){
		
			$case_id[$row_counter] = $lines[5];
			
			$patient_methylation_file_association{$case_id[$row_counter]} = $lines[1];

			$row_counter++;
			
		}
		
		@lines =();
	}

}

close INPUT;

# initialize an empty set FULL_SET_OF_GENE_NAMES
my $FULL_SET_OF_GENE_NAMES = Set::Object->new();
my @genes_per_probe;
my $test_counter = 0;
my @gene_symbol;

# for each id, get name of methylation file	
foreach my $key (keys %patient_methylation_file_association){
	# open methylation file
	my $counter=0;
	open DATA, $patient_methylation_file_association{$key} or die "Can not find file";
	#print "\n ----for $patient_methylation_file_association{$key}----\n";
	
	while (my $row = <DATA>){
		
		if ($. > 1){
			
			my @lines = split ("\t", $row);
			$genes_per_probe[$counter] = $lines[5];
			
			@gene_symbol = split (/;/, $genes_per_probe[$counter]);
			#print "hey 1 @gene_symbol \n";
			$counter++;
			$FULL_SET_OF_GENE_NAMES -> insert (@gene_symbol);
		}
	
	}
	
	
	#print "@gene_symbol/n";
	# update FULL_SET_OF_GENE_NAMES with the new name
	#print "\n hey 2 @gene_symbol\n";
	#$FULL_SET_OF_GENE_NAMES -> insert (@gene_symbol);
	close DATA;

	@genes_per_probe = ();
	@gene_symbol = ();
}

# -- at this point we have the full list of genes
#print "$FULL_SET_OF_GENE_NAMES\n";

# init the case list METHYLATION_INFO_OF_CASES
# for each id
	# get name of methylation file
	# exractMethylationVector(methylationFileName, FULL_SET_OF_GENE_NAMES) # Returns the methylation vector of the case with given ID, using the full list of genes extracted above
		# the methylation vector is a mapping between gene and methylation beta value
	# add returned vector to METHYLATION_INFO_OF_CASES
# store the METHYLATION_INFO_OF_CASES to a file
##################
# for each id, get name of methylation file	

my $test=0;
foreach my $key (keys %patient_methylation_file_association){
	
	exractMethylationVector($patient_methylation_file_association{$key},$FULL_SET_OF_GENE_NAMES);

}




# function exractMethylationVector(filename, fullListOfGeneNames)
	# create a map RES_VEC containing all the gene names as keys and NA as values
	# open file
	# for every gene name I find in the file
		# read the corresponding value from the file
		# update the RES_VEC vector for the specific gene with the value we just read
	# return RES_VEC

	
sub exractMethylationVector{

	my ($filename, $fullLSetOfGeneNames) = @_;
	my $size = $fullLSetOfGeneNames -> size();
	my @set_elements =  @$fullLSetOfGeneNames;
	my @beta_value;
	my @column;
	my @genes_probe;
	my @gene;
	my $counter;
	my $beta_counter = 0;
	print "exractMethylationVector is called!\n";
	# create a hash RES_VEC containing all the gene names as keys and NA as values
	my %RES_VEC;
	for ($counter=0; $counter < $size; $counter++){
		
		$RES_VEC{$set_elements[$counter]} = undef;
	
	}
	
	#open file
	
	print "\n$filename\n";
#	print "\n size $size\n";
#sygkrisi olwn twn genes pou exoun apothikeutei sto set me oles tis grammes gia ola ta arxeia
#sygkrisi kathe gene sto @set_elements me oles tis grammes kathe arxeioy
	for ( $counter=0; $counter < $size; $counter++){
		
		open FILE, $filename or die "Cannot find file";
		#print "\n heyyyy $set_elements[$counter]\n";
		while (my $row = <FILE>){
			if ($. > 1){
			#print "geiaa\n";
				 @column = split ("\t", $row);
			#	print "\n$row \n";
				$genes_probe[$counter] = $column[5];
				
				 @gene = split (/;/, $genes_probe[$counter]);
				#print "hey 1 @gene_symbol \n";
				#print "\n@ gene @gene\n";
				my $length_gene = @gene;
			#	print "there you go $length_gene\n";
				#print "OK! $beta_value[$] equals $column[1] \n";
				
				for (my $gene_array =0; $gene_array < $length_gene; $gene_array++){
			#	print "\nCompare $set_elements[$counter] with $gene[$gene_array]\n";	
					if ($set_elements[$counter] eq $gene[$gene_array]){
						#if ($column[1] != 'N/A'){
						
							$beta_value[$beta_counter] =  $column[1]; #o pinakas apothikeuei ola ta beta value gia to idio gonidio an uparxei se parapanw apo mia seira
							$beta_counter++;
							print "OK! $set_elements[$counter] equals $gene[$gene_array] and beta_value is $beta_value[$gene_array] \n";
							$gene_array = $length_gene; #yparxei i periptwsi to $genes per probe na einai tis morfis gene1;gene1;gene1, epomenws o pinakas @gene tha ginei gene1 gene1 gene1 kai etsi to average tha einai lathos
							#print "\n$beta_value[$gene_array]\n";
							
						#}
						
					}
				
				}
			
				
			}
			
		}
		print " \n this gene of geneset $set_elements[$counter] has these @beta_value beta values\n";
		$RES_VEC{$set_elements[$counter]} = average_beta_value_per_gene(\@beta_value);
		print "\n$set_elements[$counter] has average $RES_VEC{$set_elements[$counter]}\n";
		@column = ();
		@genes_probe = (); 
		@gene = ();
		@beta_value = ();
		$beta_counter = 0;
		close FILE;
	}
	#return %RES_VEC;
	
	}
#osa exoun NA beta_value moy bgazei pws exei 0 average. to opoio einai lathos! na brw tropo na to diorthwsw!
#den thelw oi NA values na upologzontai sto average giati to 0 exei fysiki simasia
 sub average_beta_value_per_gene{
	
	no warnings;
	
	my ($ref_beta_value) = @_;
	my @beta_value_def = @{$ref_beta_value};
	return unless @beta_value_def;
	my $length = @beta_value_def;
	my $total;
	my @beta_value_array;
	for (my $counter =0; $counter < $length; $counter ++){
		
		if ($beta_value_def[$counter] != 'N/A'){
		
			$beta_value_array[$counter] = $beta_value_def[$counter];
		}
		
	}
	foreach (@beta_value_array) {
	
			$total += $_;
		
        
    }
    return $total/scalar @beta_value_array;
}