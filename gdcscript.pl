#!usr/bin/perl -w

use strict;
use warnings;
use 5.010;
use Data::Dumper;
use Set::Object;
use Set::Object qw(set);

# open mapping file
open INPUT, "gdc_sample_sheet.2018-06-06.tsv" or die "Cannot find file";



my @sample_id=();
my @folder_id = ();
my @file_id = ();
my @data_type = (); #apothikeuei ola ta diaforetika eidi dedomenwn gia to kathe deigma
my $row_counter=0; 
my @path;
# initialize an empty hash patient_file_association, value is case id and key is the associated folder and file which is unique 
#works ok!
my %patient_meth_file_association;
my %patient_mirna_file_association;
my %patient_rnaseq_file_association;
my %patient_data_association;
my %patient_type_of_data_association;
my $FULL_SET_OF_DATA_TYPES = Set::Object->new();
my $FULL_SET_OF_SAMPLES = Set::Object->new();

while (my $line = <INPUT>){
	
	if ($. > 1){
	
		my @column = split ("\t", $line);
	
			
		$folder_id[$row_counter] = $column[0];
		
		$file_id[$row_counter] = $column[1];
		
		$data_type[$row_counter] = $column[3];
		
		$sample_id[$row_counter] = $column[6];
		
		
		$patient_type_of_data_association {$sample_id[$row_counter]} = undef;
		
		$row_counter++;

		
		
		
		@column =();
	}
	
}

$FULL_SET_OF_DATA_TYPES -> insert (@data_type);
$FULL_SET_OF_SAMPLES -> insert (@sample_id);

close INPUT;


my @SET_OF_Sample_IDs =  @$FULL_SET_OF_SAMPLES;
my @datatype = ();
	
foreach my $sample (@SET_OF_Sample_IDs){
	
	open INPUT, "gdc_sample_sheet.2018-06-06.tsv" or die "Cannot find file";
	while (my $row = <INPUT>){

		my @column = split ("\t", $row);

		my $data_type = $column[3];
		my $sample_id = $column[6];
		if ($sample eq $sample_id){
		
				push (@datatype, $data_type);
		
		}

	
	}
	my $datatypes_conc = join (";", @datatype);
	$patient_type_of_data_association {$sample} = $datatypes_conc;
	close INPUT;
	@datatype = ();
}

foreach my $sample (keys %patient_type_of_data_association){
	
	my $meth =0;
	my $rnaseq = 0;
	my $mirna =0;

	open INPUT, "gdc_sample_sheet.2018-06-06.tsv" or die "Cannot find file";
	while (my $line = <INPUT>){
	
		if ($. > 1){
			
			my @column = split ("\t", $line);
			my $folder_id = $column[0];
			my $file_id = $column[1];
			my $sample_id = $column[6];
			my $data_type = $column[3];
			
			if ($sample_id eq $sample){
				
				if ( ($data_type eq 'Methylation Beta Value')&& $meth<=1) {
					
					my $path  = join ("//", $folder_id,$file_id);
					
					$patient_meth_file_association{$sample_id} = $path;
				
					$meth++;
				}

				elsif (($data_type eq 'miRNA Expression Quantification')&& $mirna<=1){
					
					my $path  = join ("//", $folder_id,$file_id);
					
					$patient_mirna_file_association{$sample_id} = $path;
					
					$mirna++;
				}

				elsif (($data_type eq 'Gene Expression Quantification')&& $rnaseq<=1){
					
					my $path  = join ("//", $folder_id,$file_id);
					
					$patient_rnaseq_file_association{$sample_id} = $path;
					
					$rnaseq++;
					
				}
					
			}
			
			@column =();
		}
	
	}
	
	close INPUT;
	$meth =0;
	$mirna =0;
	$rnaseq =0;
}


print "--------rna----\n";
print Dumper \%patient_rnaseq_file_association;
print "--------meth----\n";
print Dumper \%patient_meth_file_association;
print "--------mirna----\n";
print Dumper \%patient_mirna_file_association;


get_all_data_for_all_patients(\%patient_rnaseq_file_association,\%patient_meth_file_association,\%patient_mirna_file_association,$FULL_SET_OF_DATA_TYPES,$FULL_SET_OF_SAMPLES);

#input ena hash poy tha sysxetizei astheneis me type of data

sub get_all_data_for_all_patients {

	no warnings;
	
	my($patient_mrna_file_association_ref,$patient_meth_file_association_ref,$patient_mirna_file_association_ref, $FULL_SET_OF_DATA_TYPES,$FULL_SET_OF_SAMPLES) = @_;
	my @set_of_data_types = @$FULL_SET_OF_DATA_TYPES;
	my @set_of_samples = @$FULL_SET_OF_SAMPLES;
	my %patient_mrna_file_association = %$patient_mrna_file_association_ref;
	my %patient_meth_file_association = %$patient_meth_file_association_ref;
	my %patient_mirna_file_association = %$patient_mirna_file_association_ref;
	my $patient_meth_data_association_hash_ref;
	my $patient_mirna_data_association_hash_ref;
	my $patient_rnaseq_data_association_hash_ref;
	my $data_type_length = @set_of_data_types;
	my $gene_sorted_list;
	my $mirna_sorted_list;
	my $gene_meth_sorted_list;
	
	for(my $counter =0; $counter<$data_type_length; $counter++){

			if ($set_of_data_types[$counter] eq 'Methylation Beta Value'){
				
			($patient_meth_data_association_hash_ref, $gene_meth_sorted_list) = get_gene_symbols_from_methylation_files (\%patient_meth_file_association);
			}

			elsif ($set_of_data_types[$counter] eq 'miRNA Expression Quantification'){

				($patient_mirna_data_association_hash_ref, $mirna_sorted_list) = get_mirna_id_from_mirna_files(\%patient_mirna_file_association);
			}

			elsif ($set_of_data_types[$counter] eq 'Gene Expression Quantification'){

				($patient_rnaseq_data_association_hash_ref, $gene_sorted_list)	= get_gene_symbols_from_rnaseq_files(\%patient_mrna_file_association);
				
			}
	}
	
	my @gene_symbol_sorted = @$gene_sorted_list;
	my %patient_meth_data_association = %$patient_meth_data_association_hash_ref;
	my @mirna_sorted = @$mirna_sorted_list;
	my %patient_mirna_data_association = %$patient_mirna_data_association_hash_ref;
	my @gene_meth_sorted = @$gene_meth_sorted_list;
	my %patient_rnaseq_data_association = %$patient_rnaseq_data_association_hash_ref;
	
	print "-----meth---\n";
	print Dumper \%patient_meth_data_association;
	print "-----mirna---\n";
	print Dumper \%patient_mirna_data_association;
	print "-----mrna---\n";
	print Dumper \%patient_rnaseq_data_association;
	
	write_output(\%patient_type_of_data_association,\%patient_meth_data_association,\%patient_mirna_data_association,\%patient_rnaseq_data_association, \@gene_symbol_sorted,\@mirna_sorted,\@gene_meth_sorted,\@set_of_samples);
	
}

#sub poy tha exei ws input 4 (mazi me ta klinika) hash enos deigmatos {key gene/gene/mirna => value } kai tis listes tous
#gia kathe stoixeio tis listas bres to key kai grapse se ena arxeio to value toy hash

sub write_output{
	
	no warnings;

	my($patient_type_of_data_association_ref,$patient_meth_data_association_hash_ref,$patient_mirna_data_association_hash_ref,$patient_rna_data_association_hash_ref, $gene_sorted_list, $mirna_sorted_list, $gene_meth_sorted_list,$all_samples) = @_;
	my @gene_symbol_list= @$gene_sorted_list;
	my @mirna_symbol_list= @$mirna_sorted_list;
	my @gene_meth_symbol_list= @$gene_meth_sorted_list;
	my @set_of_all_samples= @$all_samples;
	my %patient_type_of_data_association = %$patient_type_of_data_association_ref;
	my %patient_meth_data_association = %$patient_meth_data_association_hash_ref;
	my %patient_mirna_data_association = %$patient_mirna_data_association_hash_ref;
	my %patient_rna_data_association = %$patient_rna_data_association_hash_ref;

	my $temp = 0;
	my $counter;
	
	open WRITE, ">>output.txt" or die "cannot open file";
	#print headers
	
	print WRITE "DATA\t";
	
	foreach my $element (@gene_symbol_list){

		print WRITE "m_$element\t";

	}
	foreach my $mirna (@mirna_symbol_list){

		print WRITE "$mirna\t";

	}
	foreach my $dna_meth (@gene_meth_symbol_list){

		print WRITE "$dna_meth\t";

	}
	#print WRITE "\n";
	#close WRITE;
	
	#open WRITE, ">>output.txt" or die "cannot open file";
	
	#my $header_counter = 0;
	my $gene_length = @gene_symbol_list;
	my $mirna_length = @mirna_symbol_list;
	my $meth_length = @gene_meth_symbol_list;
	#print "length $length\n";
	
	#foreach my $sample (@set_of_all_samples){
	foreach my $sample (keys %patient_type_of_data_association){
			
			print "\n patient $sample\n"; 
			my @data_type_per_sample = split(";",$patient_type_of_data_association{$sample});
			my $length = @data_type_per_sample;
			my $rnaseq = 0;
			my $mirna = 0;
			my $meth = 0;
			#my $written_once = 0;
			print WRITE "\n$sample\t";
			
			for (my $counter_type=0; $counter_type<$length; $counter_type++){
				
				if ($data_type_per_sample[$counter_type] eq 'Gene Expression Quantification'){
				
					$rnaseq++;

				}
				elsif ($data_type_per_sample[$counter_type] eq 'miRNA Expression Quantification'){
				
					$mirna++;

				}
				elsif($data_type_per_sample[$counter_type] eq 'Methylation Beta Value'){
				
					$meth++;
				}
			}	
			
			if ($rnaseq == 1){
			
				foreach my $rna_patient (keys %patient_rna_data_association){
				
					if ($rna_patient eq $sample){
						
						for ($counter = 0; $counter < $gene_length; $counter++ ){
							
							foreach my $gene ( keys %{$patient_rna_data_association{$rna_patient}}){

									if($gene_symbol_list[$counter] eq $gene){

										#print "\n $gene_symbol_list[$counter] equals $gene and has this value $patient_rna_data_association{$rna_patient}{$gene}\n";
										
										print WRITE "$patient_rna_data_association{$rna_patient}{$gene}\t";
										#$written_once++;
									}
								
							}

						}
					}
					
				
				}
			}
			elsif($rnaseq == 0){

				for ($counter = 0; $counter < $gene_length; $counter++ ){

							print WRITE "-\t";

				}
			}

			
		if($mirna == 1){
		
				foreach my $mirna_patient (keys %patient_mirna_data_association){
				
					if ($mirna_patient eq $sample){
						
						for ($counter = 0; $counter < $mirna_length; $counter++ ){
					
							foreach my $mirna ( keys %{$patient_mirna_data_association{$mirna_patient}}){
								
								#print "mirna $mirna\n";

									if($mirna_symbol_list[$counter] eq $mirna){

										#print "\n $mirna_symbol_list[$counter] equals $mirna and has this value $patient_mirna_data_association{$mirna_patient}{$mirna}\n";
										
										print WRITE "$patient_mirna_data_association{$mirna_patient}{$mirna}\t";
									}
								
							}
				
						}
					}
				}
		}			
		elsif($mirna == 0){

			for ($counter = 0; $counter < $mirna_length; $counter++ ){

					print WRITE "-\t";

			}

		}
			
		
		if($meth == 1){
		
			foreach my $meth_patient (keys %patient_meth_data_association){
			
				if ($meth_patient eq $sample){
					
					for ($counter = 0; $counter < $meth_length; $counter++ ){
				
						foreach my $meth_gene ( keys %{$patient_meth_data_association{$meth_patient}}){
							
							#print "meth_gene $meth_gene\n";

								if($gene_meth_symbol_list[$counter] eq $meth_gene){

									#print "\n $gene_meth_symbol_list[$counter] equals $meth_gene and has this value $patient_meth_data_association{$meth_patient}{$meth_gene}\n";
									
									print WRITE "$patient_meth_data_association{$meth_patient}{$meth_gene}\t";
								}
							
							
						}
			
					}
				}
		
			}	
			
		}
		elsif($mirna == 0){
			
			for ($counter = 0; $counter < $mirna_length; $counter++ ){

								print WRITE "-\t";
		
			}
		
		}

		
		
		@data_type_per_sample = ();
		

	}
	close WRITE;			
}



sub get_gene_symbols_from_rnaseq_files{
	
	no warnings;
	
	print "sub rna seq is called!\n";
	my ($hash_ref) = @_;
    my %patient_rnaseq_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_miRNAs
	my $FULL_SET_OF_ENSEMBL_IDS = Set::Object->new();
	my $FULL_SET_OF_GENE_NAME = Set::Object->new();
	my $test_counter = 0;
	my @ensembl_id = ();
	my @gene_name;
	my @gene_symbol;
	my $gene_symbol_ref;
	my $gene_symbol_m;
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
							chomp $mart_gene_name;
							chomp $ensembl_id[$counter];
							if ($ensembl_id[$counter] eq $mart_ensembl_id){
								
								$gene_symbol_m = $mart_gene_name;
								#print "Ensembl $mart_ensembl_id or $ensembl_id has gene name $mart_gene_name";
								push(@gene_name,$gene_symbol_m);
							
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


	
	#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
	my %patient_rnaseq_data_association;
	foreach my $key (keys %patient_rnaseq_file_association){
		#print "hey babe\n$patient_mirna_file_association{$key}\n";
		#my $gene_ensembl_ref;
		
		($patient_rnaseq_data_association{$key}, $gene_symbol_ref) = exractRNAseqVector($patient_rnaseq_file_association{$key},$FULL_SET_OF_ENSEMBL_IDS, $FULL_SET_OF_GENE_NAME);
		
		@gene_symbol = @$gene_symbol_ref;
	#	print "@gene_symbol/n";
		
	}
	
	#print Dumper \%patient_data_association;

	return (\%patient_rnaseq_data_association,\@gene_symbol);
	

}



sub get_mirna_id_from_mirna_files{

	no warnings;
	
	print "sub get_mirna_id_from_mirna_files is called!\n";
	my ($hash_ref) = @_;
    my %patient_mirna_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_miRNAs
	my $FULL_SET_OF_miRNAs = Set::Object->new();
	my @miRNAs_id;

	# for each case id, get name of folder and methylation file	
	foreach my $key (keys %patient_mirna_file_association){
		# open file
		my $counter=0;
		open DATA, $patient_mirna_file_association{$key} or die "Cannot find file";
		
		while (my $row = <DATA>){
			
			if ($. > 1){
				
				my @lines = split ("\t", $row);
				$miRNAs_id[$counter] = $lines[0];		
				$counter++;
				
			}
		
		}
		$FULL_SET_OF_miRNAs -> insert (@miRNAs_id);
		
		close DATA;
		#print "mirna id @miRNAs_id\n";
		@miRNAs_id = ();
}
	#print " set $FULL_SET_OF_miRNAs\n";
	#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
	my $mirna_id_sorted_ref;
	my %patient_mirna_data_association;
	my @mirna_id_sorted;
	foreach my $patient (keys %patient_mirna_file_association){
	
		($patient_mirna_data_association{$patient},$mirna_id_sorted_ref )= exractmiRNAVector($patient_mirna_file_association{$patient},$FULL_SET_OF_miRNAs);

	}
	@mirna_id_sorted = @$mirna_id_sorted_ref;
	#print Dumper \%patient_mirna_file_association;
	return (\%patient_mirna_data_association, \@mirna_id_sorted);

}



sub get_gene_symbols_from_methylation_files{
	
	no warnings;
	
	print "meth sub is called!\n";
	my ($hash_ref) = @_;
    my %patient_methylation_file_association = %$hash_ref;
		
	# initialize an empty set FULL_SET_OF_GENE_NAMES
	my $FULL_SET_OF_GENE_NAMES = Set::Object->new();
	my @genes_per_probe = ();
	my @gene_symbol = ();
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
				$counter++;
				$FULL_SET_OF_GENE_NAMES -> insert (@gene_symbol);
			}
		
		}
		
		
		
		#$FULL_SET_OF_GENE_NAMES -> insert (@gene_symbol);
		close DATA;

		@genes_per_probe = ();
		@gene_symbol = ();
	}

	# -- at this point we have the full list of genes- it works fine

	# init the case list METHYLATION_INFO_OF_CASES
	# for each id
		# get name of methylation file
		# exractMethylationVector(methylationFileName, FULL_SET_OF_GENE_NAMES) # Returns the methylation vector of the case with given ID, using the full list of genes extracted above
			# the methylation vector is a mapping between gene and methylation beta value
		# add returned vector to METHYLATION_INFO_OF_CASES
	# store the METHYLATION_INFO_OF_CASES to a file
	##################
	# for each id, get name of methylation file	
	#creation of a hash that has as a key the case id and value a hash that is returned from the exractMethylationVector subroutine
	my $gene_sorted_meth_ref;
	my %patient_meth_data_association;
	my @gene_sorted_meth;
	foreach my $patient (keys %patient_methylation_file_association){
		
		($patient_meth_data_association{$patient},$gene_sorted_meth_ref) = exractMethylationVector($patient_methylation_file_association{$patient},$FULL_SET_OF_GENE_NAMES);

	}
	@gene_sorted_meth = @$gene_sorted_meth_ref;
	print "Also Here!!!\n";
	print Dumper \%patient_meth_data_association;
	return (\%patient_meth_data_association,\@gene_sorted_meth);

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
	my $sizeEnsembl = $fullLSetOfEnsemblId -> size();
	my $sizeEntrez = $fullLSetOfEntrezGeneId -> size();
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
	
	#sort @set_ensembl and @set_gene_name alphabetically
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_ensembl;
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_gene_name;
	
	for ($counter=0; $counter < $sizeEntrez; $counter++){
		
		$RES_VEC{$set_gene_name[$counter]} = undef;
	
	}

	#sygkrisi olwn twn genes pou exoun apothikeutei sto set me oles tis grammes gia ola ta arxeia
	#sygkrisi kathe gene sto @set_elements me oles tis grammes kathe arxeioy
	for ( $counter=0; $counter < $sizeEnsembl; $counter++){
		
		open FILE, $filename or die "Cannot find file";
		#print "\n heyyyy $set_elements[$counter]\n";
		while (my $row = <FILE>){
		
		
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
							chomp $mart_gene_name;
							
							if ($ensembl_id eq $mart_ensembl_id){
								
								chomp $column[1];
								#print "Ensembl $mart_ensembl_id or $ensembl_id has gene name $mart_gene_name";
								$RES_VEC{$mart_gene_name} = $column[1];
								#chomp $RES_VEC{$mart_gene_name};
							}
						}
						
					}
					close ANN;
					
				}
			
			
		}

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
	
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_elements;
	for ($counter=0; $counter < $size; $counter++){
		
		$RES_VEC{$set_elements[$counter]} = undef;

	}


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
				$miRNA_id= $column[0];

				if ($set_elements[$counter] eq $miRNA_id){
					
						$RES_VEC{$set_elements[$counter]} = $column[2];
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
		return (\%RES_VEC,\@set_elements);
}
	
sub exractMethylationVector{

	no warnings;
	#use warnings;
	
	print "exractMethylationVector is called!\n";
	my ($filename, $fullLSetOfGeneNames) = @_;
	my $size = $fullLSetOfGeneNames -> size();
	my @set_elements =  @$fullLSetOfGeneNames;
	my @beta_value;
	my @column;
	my @genes_probe;
	my @gene;
	my $counter;
	my $beta_counter = 0;
	
	# create a hash RES_VEC containing all the gene names as keys and NA as values
	my %RES_VEC;
	
	sort { (lc($a) cmp lc($b)) or ($a cmp $b) } @set_elements;
	for ($counter=0; $counter < $size; $counter++){
		
		$RES_VEC{$set_elements[$counter]} = undef;
	
	}
	
	#open file
	
	#print "\n$filename\n";
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
	
		print "Heereee\n";
		print Dumper \%RES_VEC;
	return (\%RES_VEC,\@set_elements);	
	
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
		#print "array_length is $array_length \n";
		#for ( $counter =0; $counter< $array_length; $counter++){
		#	print "\nhey\n";
		#	$total = $beta_value_array[$counter]+ $total; 
			
		#}
		foreach (@beta_value_array) {
		
				$total += $_;
			
			
	    }

		
		return $total/$array_length;
	
	}
	
}