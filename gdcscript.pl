#!usr/bin/perl -w

use strict;
use warnings;
use 5.010;
use Data::Dumper;
use Set::Object;
use Set::Object qw(set);

# open mapping file
open INPUT, "gdc_sample_sheet.2018-06-06.tsv" or die "Cannot find file";

#print "Hello! What type of data you want to process?\n";

my $data_type = "Gene Expression Quantification";
chomp $data_type;

my @sample_id=();
my @folder_id = ();
my @file_id = ();
my @data_type = (); #apothikeuei ola ta diaforetika eidi dedomenwn gia to kathe deigma
my $row_counter=0; 
my @path;
# initialize an empty hash patient_file_association, value is case id and key is the associated folder and file which is unique 
#works ok!
my %patient_file_association;
my %patient_data_association;

while (my $line = <INPUT>){
	
	if ($. > 1){
	
		my @column = split ("\t", $line);
		
		if($column[3] eq $data_type){
			
			#print "OK\n";
			
			$folder_id[$row_counter] = $column[0];
			
			$file_id[$row_counter] = $column[1];
			
			$data_type[$row_counter] = $column[3];
			
			$sample_id[$row_counter] = $column[6];
			
			$path[$row_counter]  = join ("//", $folder_id[$row_counter],$file_id[$row_counter]);
			
			$patient_file_association{$sample_id[$row_counter]} = $path[$row_counter];
			
			$patient_data_association{$sample_id[$row_counter]} = undef;
			
			$row_counter++;
			
		}
		
		#get_all_data_for_all_patients(%patient_file_association);
		
		@column =();
	}
	#print Dumper \%patient_file_association;
}
close INPUT;

#input ena hash poy tha sysxetizei astheneis me type of data
#sub get_all_data_for_all_patients {

#	my $data_type_length= @data_type;

	#for(my $counter =0; $counter<$data_type_length; $counter++){

		#	if ($data_type[$counter] eq 'Methylation Beta Value'){
				
			#	get_gene_symbols_from_methylation_files (\%patient_file_association);

			#}

		#	elsif ($data_type[$counter] eq 'miRNA Expression Quantification'){

		#		get_mirna_id_from_mirna_files(\%patient_file_association);
		#	}

		#	elsif ($data_type[$counter] eq 'Gene Expression Quantification'){

	my ($patient_data_association_hash_ref, $gene_sorted_list)	= get_gene_symbols_from_rnaseq_files(\%patient_file_association);
	
	my @gene_symbol_sorted = @$gene_sorted_list;
	my %patient_data_association_hash = %$patient_data_association_hash_ref;
	
	print "\n Edwwww\n";
	print Dumper \%patient_data_association_hash;
	
	#write_output($patient_data_association_hash_ref, \@gene_symbol_sorted);
		#	}

#}
#}

#sub poy tha exei ws input 4 (mazi me ta klinika) hash enos deigmatos {key gene/gene/mirna => value } kai tis listes tous
#gia kathe stoixeio tis listas bres to key kai grapse se ena arxeio to value toy hash

sub write_output{
	
	my($patient_data_association_hash_ref, $gene_sorted_list) = @_;
	my @gene_symbol_list= @$gene_sorted_list;
	my %patient_data_association = %$patient_data_association_hash_ref;
	my $temp = 0;
	open WRITE, (">>output.txt" or die "cannot open file");
	#print headers
	if ($temp == 0 ){
	
		print WRITE "DATA";
		
		foreach my $element (@gene_symbol_list){
	
			print WRITE "$element";
	
		}
		$temp++;
		close WRITE;
	}
	open WRITE, (">>output.txt" or die "cannot open file");
	
	my $header_counter = 0;
	my $length = @gene_symbol_list;
	
		
	foreach my $patient (keys %patient_data_association){
		
		if ($header_counter == 0){
				
				print WRITE "$patient";
				$header_counter ++;
				
			}
		foreach my $gene (keys %{$patient_data_association{$patient}}){
		
			for ( my $counter = 0; $counter< $length; $counter++ ){
			
				if($gene_symbol_list[$counter] eq $gene){
					
					print WRITE "$patient_data_association{$patient}{$gene}";
				
				}
		
			}
		
		}

		
		$header_counter = 0;
		
	}
}
	
	#foreach my $key (keys %patient_data_association){
	
		
	
	#}



sub get_gene_symbols_from_rnaseq_files{
	print "sub rna seq is called!\n";
	my ($hash_ref) = @_;
    my %patient_rnaseq_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_miRNAs
	my $FULL_SET_OF_ENSEMBL_IDS = Set::Object->new();
	my $FULL_SET_OF_GENE_NAME = Set::Object->new();
	my $test_counter = 0;
	my @ensembl_id = ();
	my @gene_name;
	my @gene_ensembl_ref;
	# for each case id, get name of folder and methylation file	
	foreach my $key (keys %patient_rnaseq_file_association){
		# open file
		my $counter=0;

		open DATA, $patient_rnaseq_file_association{$key} or die "Cannot find file";
		
		#print "\n ----for $patient_rnaseq_file_association{$key}----\n";
		
		while (my $row = <DATA>){
				
			my @lines = split ("\t", $row);
			$ensembl_id[$counter] = $lines[0];	
			
			open ANN, "mart_export_versions.txt" or die "Cannot find file";
						
					while (my $line = <ANN>){
						
						my $mart_ensembl_id;
						my $mart_gene_name;
						
						if ($.>1){
							
							my @rows = split ("\t", $line);
							$mart_ensembl_id = $rows[0];
							$mart_gene_name = $rows[1];
							
							if ($ensembl_id[$counter] eq $mart_ensembl_id){
								
								#print "Ensembl $mart_ensembl_id or $ensembl_id has gene name $mart_gene_name";
								push(@gene_name,$mart_gene_name);
							
							}
						}
						
					}
					close ANN;
					
			#print "ensemble ids $ensembl_id[$counter]\n";			
			$counter++;
	
		}
		my $length = @ensembl_id;
	#	print "$length";
	
		$FULL_SET_OF_ENSEMBL_IDS -> insert (@ensembl_id);
		$FULL_SET_OF_GENE_NAME -> insert (@gene_name);
		close DATA;
		
		
		@ensembl_id = ();
		@gene_name = ();
}
	#my $size = $FULL_SET_OF_ENSEMBL_IDS->size();
	#print " set $size\n";
	#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
	foreach my $key (keys %patient_rnaseq_file_association){
		#print "hey babe\n$patient_mirna_file_association{$key}\n";
		my $gene_ensembl_ref;
		
		($patient_data_association{$key}, $gene_ensembl_ref) = exractRNAseqVector($patient_rnaseq_file_association{$key},$FULL_SET_OF_ENSEMBL_IDS, $FULL_SET_OF_GENE_NAME);
		
		@gene_ensembl_ref = @$gene_ensembl_ref;
		print "@gene_ensembl_ref/n";
		
	}

	return (\%patient_data_association,\@gene_ensembl_ref);
	#print Dumper \%patient_data_association;

}



sub get_mirna_id_from_mirna_files{
	print "sub is called!\n";
	my ($hash_ref) = @_;
    my %patient_mirna_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_miRNAs
	my $FULL_SET_OF_miRNAs = Set::Object->new();
	my $test_counter = 0;
	my @miRNAs_id;

	# for each case id, get name of folder and methylation file	
	foreach my $key (keys %patient_mirna_file_association){
		# open file
		my $counter=0;
		open DATA, $patient_mirna_file_association{$key} or die "Cannot find file";
		#print "\n ----for $patient_mirna_file_association{$key}----\n";
		
		while (my $row = <DATA>){
			
			if ($. > 1){
				
				my @lines = split ("\t", $row);
				$miRNAs_id[$counter] = $lines[0];		
				$counter++;
				$FULL_SET_OF_miRNAs -> insert (@miRNAs_id);
			}
		
		}
		
		close DATA;
		print "mirna id @miRNAs_id\n";
		@miRNAs_id = ();
}
	print " set $FULL_SET_OF_miRNAs\n";
	#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
	foreach my $key (keys %patient_mirna_file_association){
		#print "hey babe\n$patient_mirna_file_association{$key}\n";
		$patient_data_association{$key} = exractmiRNAVector($patient_mirna_file_association{$key},$FULL_SET_OF_miRNAs);

	}

	print Dumper \%patient_data_association;

}



sub get_gene_symbols_from_methylation_files{
	print "!!!!!!!!!!! meth sub is called!\n";
	my ($hash_ref) = @_;
    my %patient_methylation_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_GENE_NAMES
	my $FULL_SET_OF_GENE_NAMES = Set::Object->new();
	my @genes_per_probe = ();
	my $test_counter = 0;
	my @gene_symbol;

	# for each case id, get name of folder and methylation file	
	foreach my $key (keys %patient_methylation_file_association){
		# open methylation file
		my $counter=0;
		open DATA, $patient_methylation_file_association{$key} or die "Cannot find file";
		#print "\n ----for $patient_methylation_file_association{$key}----\n";
		
		while (my $row = <DATA>){
			
			if ($. > 1){
				
				my @lines = split ("\t", $row);
				$genes_per_probe[$counter] = $lines[5];		
				@gene_symbol = split (/;/, $genes_per_probe[$counter]);
			#	print "hey  @gene_symbol \n";
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

# -- at this point we have the full list of genes- it works fine
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
#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
foreach my $key (keys %patient_methylation_file_association){
	#print "hey babe\n$patient_methylation_file_association{$key}\n";
	$patient_data_association{$key} = exractMethylationVector($patient_methylation_file_association{$key},$FULL_SET_OF_GENE_NAMES);

}

print Dumper \%patient_data_association;

}



# function exractMethylationVector(filename, fullListOfGeneNames)
	# create a map RES_VEC containing all the gene names as keys and NA as values
	# open file
	# for every gene name I find in the file
		# read the corresponding value from the file
		# update the RES_VEC vector for the specific gene with the value we just read
	# return RES_VEC
	
sub exractRNAseqVector{

	no warnings;
	
	my ($filename, $fullLSetOfEnsemblId, $fullLSetOfEntrezGeneId ) = @_;
	my $size = $fullLSetOfEnsemblId -> size();
	my @set_ensembl =  @$fullLSetOfEnsemblId;
	my @set_gene_name = @$fullLSetOfEntrezGeneId;
	my $FPKM;
	my @column;
	my $miRNA_id;
	my @gene;
	my $counter;
	my $reads_counter = 0;
	my $ensembl_id;
	
	print "exractRNAseqVector is called!\n";
	# create a hash RES_VEC containing all the gene names as keys and NA as values
	my %RES_VEC;
	for ($counter=0; $counter < $size; $counter++){
		
		$RES_VEC{$set_gene_name[$counter]} = undef;
	
	}
	
	#sort @set_ensembl and @set_gene_name alphabetically
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_ensembl;
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_gene_name;
	
#	for ( $counter=0; $counter < $size; $counter++){
		
	#	print "\nsort ensembl $set_ensembl[$counter]\n ";
	
	#	}
		
	#	for ( $counter=0; $counter < $size; $counter++){
		
	#	print "\nsort gene_symbol $set_gene_name[$counter]\n ";
	
		#}

	
	#open file
	
	#print "\n$filename\n";
	#	print "\n size $size\n";
	#sygkrisi olwn twn genes pou exoun apothikeutei sto set me oles tis grammes gia ola ta arxeia
	#sygkrisi kathe gene sto @set_elements me oles tis grammes kathe arxeioy
	for ( $counter=0; $counter < $size; $counter++){
		
		open FILE, $filename or die "Cannot find file";
		#print "\n heyyyy $set_elements[$counter]\n";
		while (my $row = <FILE>){
			
			if ($. > 1){
		
				@column = split ("\t", $row);
			
				$ensembl_id= $column[0];
					
				if ($set_ensembl[$counter] eq $ensembl_id){
					
					open ANN, "mart_export_versions.txt" or die "Cannot find file";
						
					while (my $line = <ANN>){
						
						my $mart_ensembl_id;
						my $mart_gene_name;
						
						if ($.>1){
							
							my @rows = split ("\t", $line);
							$mart_ensembl_id = $rows[0];
							$mart_gene_name = $rows[1];
							
							if ($ensembl_id eq $mart_ensembl_id){
								
								#print "Ensembl $mart_ensembl_id or $ensembl_id has gene name $mart_gene_name";
								$RES_VEC{$mart_gene_name} = $column[1];
							
							}
						}
						
					}
					close ANN;
					
				}
			
				}
		}
		#print " \n this gene of geneset $set_elements[$counter] has these @beta_value beta values\n";
		
		#$RES_VEC{$set_elements[$counter]} = $reads_per_million_miRNA_mapped ;
		#print "\n$set_elements[$counter] has $RES_VEC{$set_elements[$counter]}\n";
		@column = ();
		#@reads_per_million_miRNA_mapped = ();
		#$reads_counter = 0;
		close FILE;
		}
		#print Dumper \%RES_VEC;
		return (\%RES_VEC, \@set_gene_name);
	}
	
	
	sub exractmiRNAVector{

	no warnings;
	
	my ($filename, $fullLSetOfmiRNAs) = @_;
	my $size = $fullLSetOfmiRNAs -> size();
	my @set_elements =  @$fullLSetOfmiRNAs;
	my $reads_per_million_miRNA_mapped;
	my @column;
	my $miRNA_id;
	my @gene;
	my $counter;
	my $reads_counter = 0;
	print "exractmiRNAVector is called!\n";
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
				#print "\n$row \n";
				$miRNA_id= $column[0];
				

			if ($set_elements[$counter] eq $miRNA_id){
				
					#$reads_per_million_miRNA_mapped =  $column[2]; 
					$RES_VEC{$set_elements[$counter]} = $column[2];
			}
		
			}
		}
		#print " \n this gene of geneset $set_elements[$counter] has these @beta_value beta values\n";
		
		#$RES_VEC{$set_elements[$counter]} = $reads_per_million_miRNA_mapped ;
		print "\n$set_elements[$counter] has $RES_VEC{$set_elements[$counter]}\n";
		@column = ();
		#@reads_per_million_miRNA_mapped = ();
		#$reads_counter = 0;
		close FILE;
		}
		print Dumper \%RES_VEC;
		return \%RES_VEC;
	}
	
sub exractMethylationVector{

	no warnings;
	
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
							#print "OK! $set_elements[$counter] equals $gene[$gene_array] and beta_value is $beta_value[$gene_array] \n";
							$gene_array = $length_gene; #yparxei i periptwsi to $genes per probe na einai tis morfis gene1;gene1;gene1, epomenws o pinakas @gene tha ginei gene1 gene1 gene1 kai etsi to average tha einai lathos
							#print "\n$beta_value[$gene_array]\n";
							
						#}
					}
				}
			}
		}
		#print " \n this gene of geneset $set_elements[$counter] has these @beta_value beta values\n";
		
		$RES_VEC{$set_elements[$counter]} = average_beta_value_per_gene(\@beta_value);
		print "\n$set_elements[$counter] has average $RES_VEC{$set_elements[$counter]}\n";
		@column = ();
		@genes_probe = (); 
		@gene = ();
		@beta_value = ();
		$beta_counter = 0;
		close FILE;
		}
		#print Dumper \%RES_VEC;
	return \%RES_VEC;
	}
#osa exoun NA beta_value moy bgazei pws exei 0 average. to opoio einai lathos! na brw tropo na to diorthwsw!
#den thelw oi NA values na upologzontai sto average giati to 0 exei fysiki simasia
 sub average_beta_value_per_gene{
	
	no warnings;
	#print "\naverage_beta_value_per_gene is called!\n";
	my ($ref_beta_value) = @_;
	my @beta_value_def = @{$ref_beta_value};
	#print "\n beta_value_def is @beta_value_def\n ";
	return unless @beta_value_def;
	my $length = @beta_value_def;
	#print "length $length\n";
	my $total;
	my $array_length; 
	my @beta_value_array = ();

	my $counter;
	for ( $counter =0; $counter < $length; $counter ++){
		
		if ($beta_value_def[$counter] != 'N/A'){
			print "\ngeiaa\n";
			$beta_value_array[$counter] = $beta_value_def[$counter];
			$array_length++;
		}
		
	}
	#print "\n beta_value_array is @beta_value_array\n ";
	if ($array_length == 0 ){
		
		return undef;
	
	}
	
	else{
		my $total =0;
		#my $array_length = @beta_value_array;
		print "array_length is $array_length \n";
		#for ( $counter =0; $counter< $array_length; $counter++){
		#	print "\nhey\n";
		#	$total = $beta_value_array[$counter]+ $total; 
			
		#}
		foreach (@beta_value_array) {
		
				$total += $_;
			
			
	   }

		#print "\n total $total\n";
		#print "\n array_length $array_length\n";
		
		return $total/$array_length;
	
	}
	
}