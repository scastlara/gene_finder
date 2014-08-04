#!/usr/bin/perl
#
# Article_downloader.pl - This script reads a website that contains a 
#						  list of gene symbols and its corresponding 
#						  synonyms and it creates a synonyms hash table
#
#			Arguments: 	  none
#
# 			Requirements: none
#
#********************************************************************************


use warnings;
use strict;
use Data::Dumper;
use LWP::Simple;


#********************************************************************************
# VARIABLES 
#********************************************************************************

my $text = "";
my %hash_table = ();
my $synonyms = shift @ARGV;
#my $problem_file = shift @ARGV;




#********************************************************************************
# MAIN LOOP 
#********************************************************************************


# Creating hash table...

&table_reader($synonyms);

print Data::Dumper-> Dump ([ \ %hash_table,], [ qw/ *HASH TABLE/ ]) ;


# Looking for genes in text...

#open (TEXTFILE, $file);

#while (my $line = <TEXTFILE>) {

#	chomp $line;
	
#	my @words = split /\s/, $line;
	
#	@words = map { $_=~ s/[^A-Z0-9\s]//gi;
#				   uc $_ } @words;

#	foreach my $word (@words) {


#		if (exists $hash_table{$word}) {


#			if ($hash_table{$word}->[0] == 0) {

#				print "$word : $hash_table{$word}->[2]\n";
#
#			} # if index = 0





#		} # if primary key exists





#	} # foreach word
	




#} # while <TEXTFILE>







#********************************************************************************
# FUNCTIONS 
#********************************************************************************

#********************************************************************************
# table_reader()
#
# Arguments: webpage with symbols and synonyms in plain text
#
# Returns: nothing, it creates an array with synonyms.
#		   it calls hash_creator()
#
#


sub table_reader($) {
	my $synonyms = shift;

	# reading webpage line by line

	open (SYN, $synonyms);
	
	<SYN>; # skip 1st line
	
	while (my $line = <SYN>) {
		
		next if ($line =~ m/withdrawn/g);
		my @fields = split /\t/, $line;
		my @all_synonyms = ();
		my @normal_fields = ();
		my @quoted_fields = ();
		
		my $Ap_SYMBOL = $fields[0];
		
		if ($fields[1]) {

			my $Ap_NAME = $fields[1];
			push @normal_fields, $Ap_NAME;

		}; 

		if ($fields[2]) {

			my @Prev_SYMBOL = split /,/, $fields[2];
			@Prev_SYMBOL = map {$_=~ s/^\s//; $_} @Prev_SYMBOL;
			push @normal_fields, @Prev_SYMBOL;
		}; 


		if ($fields[3]) {

			my @SYN = split /,/, $fields[3];
			@SYN = map {$_=~ s/^\s//; $_} @SYN;
			push @normal_fields, @SYN;
		};

		
		if ($fields[4]) {

			my $ENSEMBL = $fields[4];
			push @normal_fields, $ENSEMBL;
		}; 

		
		if ($fields[5]) {

			my $REF_SEQ = $fields[5];
			push @normal_fields, $REF_SEQ;
		}; 

		if ($fields[6]) {

			my $UniProt = $fields[6];
			push @normal_fields, $UniProt;

		};

		if ($fields[7]) {

			my $prev_NAMES = $fields[7];

			while ($prev_NAMES =~ m/"(.+?)"/g) {

				my $match = $1;

				$match =~ s/^\s//;
    
    			push @quoted_fields, $match;
    
   			 }; # while match

		};

		if ($fields[8]) {

			my $NAME_syn = $fields[8];

			while ($NAME_syn =~ m/"(.+?)"/g) {
    
    			my $match = $1;

				$match =~ s/^\s//; # remove first character if it is a whitespace
    			push @quoted_fields, $match;
    
   			 }; # while match

		};


		@all_synonyms = (@normal_fields, @quoted_fields);


		&hash_creator($Ap_SYMBOL, \@all_synonyms);		

	} # while <SYN>

	close (SYN);


} # sub table_reader



#********************************************************************************
# hash_creator()
#
# Arguments: gene symbol
#			 ref to array with all synonyms of that symbol			 
#
# Returns: nothing, it creates a hash with the synonyms and the corresponding symbols
#          it calls words_ladder() if necessary
#
#


sub hash_creator {
	
	my $SYMBOL = shift;
	my $synonyms_ary = shift;


	foreach my $synonym (@$synonyms_ary) {

		my $ucsynoym = uc $synonym; # make synonyms uppercase
		
		$ucsynoym =~ s/[^A-Z0-9\s]//g; # remove non word/number/space characters
		
		my @syn_words = split /\s/, $ucsynoym; # get each word of synonym
		
		my $first_word = shift @syn_words;
		
		return if (!$first_word);


		if (@syn_words == 0) {

			
			$hash_table{$first_word} = [] if (!exists $hash_table{$first_word}); # initialize if doesn't exist
			
			
			if (!defined $hash_table{$first_word}->[0]) {

				$hash_table{$first_word}->[0] = '0'; # index = 0 if it's a dead end and synonym didn't exist before
			
			} elsif ($hash_table{$first_word}->[0]) {

				$hash_table{$first_word}->[0] = '2'; # index = 2 if it's a dead end BUT the index was 1

			}; # if
			
			
			$hash_table{$first_word}->[1] = undef if (!defined $hash_table{$first_word}->[1]); # undefs 2nd value (hash) if it didn't exist before
			
			$hash_table{$first_word}->[2] = $SYMBOL; # puts symbol if it's a dead end

			$hash_table{$first_word}->[3] = $synonym; # puts original synonym if it's a dead end
			
 			

		} else {

			
			$hash_table{$first_word} = [] if (!exists $hash_table{$first_word});
			
			
			&words_ladder($SYMBOL, \@syn_words, $hash_table{$first_word}, $synonym);
			
 			
	    } # if 



	} # foreach


} # sub hash_creator



#********************************************************************************
# words_ladder()
#
# Arguments: gene symbol
#			 ref to array with all the words of a synonym
#            empty hash			 
#
# Returns:   nothing
#			 it keeps calling itself until it runs out of words
#
#

sub words_ladder {
	

	my $SYMBOL = shift;
	my $syn_words_ary = shift;
	my $hash_table_hsh = shift;
	my $synonym = shift;

	my $next_word = shift @$syn_words_ary;
	
	
	if ($hash_table_hsh->[2]) {

		$hash_table_hsh->[0] = '2'; # index = 2 if word already has a symbol

	} else {
		
		$hash_table_hsh->[0] = '1'; # otherwise, index = 1

	}; # if symbol defined
	
	$hash_table_hsh->[1] = {} unless (defined $hash_table_hsh->[1]); # creates hash at position 2 if it's not defined
	$hash_table_hsh->[2] = undef unless ($hash_table_hsh->[2]); # undefs symbol if it doesn't exist (not true)
	$hash_table_hsh->[3] = undef unless ($hash_table_hsh->[3]); # undefs original synonym if it doesn't exist (not true)

	$hash_table_hsh->[1]->{$next_word} = [] if (!exists $hash_table_hsh->[1]->{$next_word});

	

	my $new_hsh = $hash_table_hsh->[1]->{$next_word};

	
	if (@$syn_words_ary == 0 and defined $new_hsh->[1]) {

		$new_hsh->[0] = '2'; # index = 2 if it's a dead end but there's a defined hash in [1]
		$new_hsh->[2] = $SYMBOL;
		$new_hsh->[3] = $synonym; # puts original synonym if it's a dead end
		return;


	} elsif (@$syn_words_ary == 0) {

		$new_hsh->[0] = '0'; # otherwise, index = 0
		$new_hsh->[1] = undef;
		$new_hsh->[2] = $SYMBOL;
		$new_hsh->[3] = $synonym; # puts original synonym if it's a dead end
		return;

	}; # if 

	&words_ladder($SYMBOL, $syn_words_ary, $new_hsh, $synonym);

}; # sub words_ladder