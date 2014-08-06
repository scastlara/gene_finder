#!/usr/bin/perl
#
# Article_downloader.pl - This script reads a file that contains a 
#						  list of gene symbols and its corresponding 
#						  synonyms and it creates a synonyms hash table.
#						  Then, it looks for matches in a text file and 
#						  saves the lines containing them as matches_originalfile.txt
#
#			Arguments: 	  Synonyms tabular file
#						  Text file to analyze
#						  Stopwords file
#
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

die "\nYou have to introduce 3 files as command line arguments:\n" .
	"\t- Gene synonyms table\n" .
	"\t- The text file you want to analyze\n" .
	"\t- Stopwords file\n\n"

	unless (@ARGV == 3);

my $synonyms = shift @ARGV;
my $problem_file = shift @ARGV;
my $stop_words_file = shift @ARGV;
my %stop_words = ();




#********************************************************************************
# MAIN LOOP 
#********************************************************************************


print STDERR "\n## STARTING PROGRAM ##\n## CREATING SYNONYMS HASH TABLE...\n";

&table_reader($synonyms);

#print STDERR Data::Dumper-> Dump ([ \ %hash_table,], [ qw/ *HASH TABLE/ ]) ;


print STDERR "\n## READING STOPWORDS LIST...\n";

&stop_words_reader($stop_words_file, \%stop_words);


# Creating output file name

my $matches_out = "matches_$problem_file"; 


print STDERR "\n## LOOKING FOR GENES IN TEXT...\n";

&text_analyzer($problem_file, \%hash_table, \%stop_words, $matches_out);


print STDERR "\n## MATCHES SAVED AS $matches_out.\n## PROGRAM FINISHED ##\n";




#********************************************************************************
# FUNCTIONS 
#********************************************************************************

#********************************************************************************
# table_reader()
#
# Arguments: file with symbols and synonyms in plain text (tbl)
#
# Returns: nothing, it creates an array with synonyms
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

	
	# SYMBOL KEYS!
	
	if (!exists $hash_table{$SYMBOL}) {

		$hash_table{$SYMBOL} = [];
		$hash_table{$SYMBOL}->[0] = '0';
		$hash_table{$SYMBOL}->[1] = undef;
		$hash_table{$SYMBOL}->[2] = $SYMBOL;
		$hash_table{$SYMBOL}->[3] = $SYMBOL; 


	} elsif ($hash_table{$SYMBOL}->[0] == 1){

		$hash_table{$SYMBOL}->[0] = '2';
		$hash_table{$SYMBOL}->[2] = $SYMBOL;
		$hash_table{$SYMBOL}->[3] = $SYMBOL; 

	} # if 



	# SYNONYM KEYS!

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
#			 it keeps calling itself until syn_words_ary runs out of words
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




#********************************************************************************
# stop_words_reader()
#
# Arguments: file with english stopwords
#			 reference to stopwords hash
#			 
# Returns:   nothing
#			 it creates stopwords keys 
#
#


sub stop_words_reader {
	

	my $file = shift;
	my $stop_words_hsh = shift;

	open (WORDS, $file) or die "Can't open stopwords file\n";
	
	while (<WORDS>) {
		
		chomp;

		next if (/^\s+/g); # skip blank lines

		$stop_words_hsh->{$_} = undef if (!exists $stop_words_hsh->{$_});


	}; # while <WORDS>


}; # sub stop_words_reader




#********************************************************************************
# text_analyzer()
#
# Arguments: file to analyze
#			 synonyms hash table
#			 stopwords hash
#			 output file name
#			 
# Returns:   nothing
#			 it finds genes in a text file and saves lines that contain them
#			 it calls recursive_search()
#
#

sub text_analyzer($$) {
	
	my $file = shift;
	my $hash = shift;
	my $stop_words_hsh = shift;
	my $outfile = shift;

	open (TEXTFILE, $file);

	open (OUTFILE, "> $outfile");


	while (my $line = <TEXTFILE>) {

		chomp $line;
	
		my @allwords = split /\s/, $line; 
	
		@allwords = map { $_=~ s/[\.,]//g;
						  $_ } @allwords; # removes commas and periods for processing

		my @words = grep { !exists $stop_words_hsh->{$_} } @allwords; # Removes stopwords for processing

		@words = map { $_=~ s/[^A-Z0-9\s]//gi;
									    uc $_ } @words; # removes not word/number characters and makes everything uppercase
		
		

		my @words_copy = @words;

		for (my $i = 0; $i < @words;) {
			
		
			my $complete_gene = "";

			&recurive_search($words[$i], \$complete_gene, $hash, $line, *OUTFILE, @words_copy);
			
			
			# skip complete gene from iteration (to avoid internal matches)
			
			if ($complete_gene) {

				my $count = $complete_gene =~ s/((^|\s)\S)/$1/g;

				$i += $count; 
			
				splice @words_copy, 0, $count;

			} else {

				$i++;
				splice @words_copy, 0, 1;

			} # if

		} # foreach word




	} # while <TEXTFILE>


	close (TEXTFILE);
	close (OUTFILE);


}; # sub text_analyzer



#********************************************************************************
# recursive_search()
#
# Arguments: first word to analyze
#			 ref to scalar with gene name (it becomes longer as the function calls itself)
#			 hash with synonyms
#			 line that it is processing
#			 output filehandle
#			 a copy of all the words in the sentence (from the first word onwards)
# 
#			 
# Returns:   nothing 
#			 it prints gene symbol, synonym and the line that contains it
#
#

sub recurive_search {
	
	my $word = shift;
	my $complete_gene = shift;
	my $hash = shift;
	my $line = shift;
	my $OUTFILE = shift;
	my @words_copy = @_;
	

	return unless ($word); # ends function if it runs out of words

	

	if (exists $hash->{$word}) {


		if ($hash->{$word}->[0] == 0) {

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} # if gene name has one word or else
			
			print $OUTFILE "MATCH ($$complete_gene : $hash->{$word}->[2]) at:\n".
				  "$line\n\n";	
				  		
			return;
	
		} elsif ($hash->{$word}->[0] == 1) {

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} # if gene name has one word or else


			splice @words_copy, 0, 1;
			
			my $new_hash = $hash->{$word}->[1];
			
			&recurive_search($words_copy[0], $complete_gene, $new_hash, $line, $OUTFILE, @words_copy);


		} elsif ($hash->{$word}->[0] == 2) {

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} # if gene name has one word or else

			

			my $new_hash = $hash->{$word}->[1];
			splice @words_copy, 0, 1;

			if ($words_copy[0] and exists $new_hash->{$words_copy[0]}) {

				&recurive_search($words_copy[0], $complete_gene, $new_hash, $line, $OUTFILE, @words_copy);

			} else {

				print $OUTFILE "MATCH ($$complete_gene : $hash->{$word}->[2]) at:\n".
				      "$line\n\n";
			
				return;

			} # if next word exists in hash table


		}; # if index = (0 or 1 or 2)



	return;

	} # if primary key exists


}; # sub recurive_search