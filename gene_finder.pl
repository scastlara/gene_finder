#!/usr/bin/perl
#
# 		 gene_finder.pl - This script reads a file that contains a 
#						  list of gene symbols and its corresponding 
#						  synonyms and it creates a synonyms hash table.
#						  Then, it looks for matches in a text file and 
#						  saves the lines containing them as matches_originalfile.txt
#
#			Arguments: 	  Synonyms tabular file
#						  Stopwords file
#						  Text file(s) to analyze
#
#
#********************************************************************************


use warnings;
use strict;
use Data::Dumper;
use LWP::Simple;

die "\nYou have to introduce at least 3 files as command line arguments:\n" .
	"\t- Gene synonyms table\n" .
	"\t- Stopwords file\n" .
	"\t- Text file(s) you want to analyze\n\n"

	unless (@ARGV >= 3);

#********************************************************************************
# VARIABLES 
#********************************************************************************

my $synonyms 		= shift @ARGV;
my $stop_words_file = shift @ARGV;
my @problem_files   = @ARGV;
my %hash_table      = ();
my %stop_words 		= ();


#********************************************************************************
# MAIN LOOP 
#********************************************************************************

print STDERR "\n## STARTING PROGRAM ##\n\n" .
			 "# CREATING SYNONYMS HASH TABLE...\n\n";

&table_reader($synonyms);

print STDERR "# READING STOPWORDS LIST...\n\n";

&stop_words_reader($stop_words_file, \%stop_words);

print STDERR "# LOOKING FOR GENES IN TEXT...\n";

foreach my $file (@problem_files) {

	print STDERR "\t# Analyzing $file...\n";

	my @tagged_lines = &text_analyzer($file, \%hash_table, \%stop_words);
	
	$file =~ s/.+\///;
	my $matches_out  = "matches_$file"; # Creating output file name
	&tagged_lines_filter(\@tagged_lines, $matches_out);

	print STDERR "\t# Matches saved as: $matches_out\n\n";

}

print STDERR "\n## PROGRAM FINISHED ##\n";


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


sub table_reader {
	
	my $synonyms = shift;

	# reading webpage line by line

	open (SYN, $synonyms)
		or die "Can't open $synonyms : $!\n";
		
	<SYN>; # skip 1st line
	
	while (my $line = <SYN>) {
		
		next if ($line =~ m/withdrawn/g);

		my @fields        = split /\t/, $line;
		my @all_synonyms  = ();
		my @normal_fields = ();
		my @quoted_fields = ();
		my $Ap_SYMBOL     = $fields[0];
		
		if ($fields[1]) {

			my $Ap_NAME = $fields[1];
			push @normal_fields, $Ap_NAME;

		}; 

		if ($fields[2]) {

			my @Prev_SYMBOL = split /,/, $fields[2];
			@Prev_SYMBOL    = map {$_=~ s/^\s//; $_} @Prev_SYMBOL;
			push @normal_fields, @Prev_SYMBOL;
		}; 

		if ($fields[3]) {

			my @SYN = split /,/, $fields[3];
			@SYN    = map {$_=~ s/^\s//; $_} @SYN;
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
				$match    =~ s/^\s//;
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

		@all_synonyms = grep { length($_) > 2} @all_synonyms;

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
#          it calls recursive_hash() if necessary
#
#


sub hash_creator {
	
	my $SYMBOL       = shift;
	my $synonyms_ary = shift;

	# SYMBOL KEYS!
	
	if (!exists $hash_table{$SYMBOL}) {

		$hash_table{$SYMBOL}      = [];
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

		my $ucsynoym   = uc $synonym; # make synonyms uppercase
		$ucsynoym      =~ s/[^A-Z0-9\s]//g; # remove non word/number/space characters
		my @syn_words  = split /\s/, $ucsynoym; # get each word of synonym
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
			&recursive_hash($SYMBOL, \@syn_words, $hash_table{$first_word}, $synonym);
	 		
	    } # if 

	} # foreach

} # sub hash_creator


#********************************************************************************
# recursive_hash()
#
# Arguments: gene symbol
#			 ref to array with all the words of a synonym
#            empty hash			 
#
# Returns:   nothing
#			 it keeps calling itself until syn_words_ary runs out of words
#
#

sub recursive_hash {
	
	my $SYMBOL         = shift;
	my $syn_words_ary  = shift;
	my $hash_table_hsh = shift;
	my $synonym        = shift;
	my $next_word      = shift @$syn_words_ary;
		
	if ($hash_table_hsh->[2]) {

		$hash_table_hsh->[0] = '2'; # index = 2 if word already has a symbol

	} else {
		
		$hash_table_hsh->[0] = '1'; # otherwise, index = 1

	}; # if symbol defined
	
	$hash_table_hsh->[1] = {} unless (defined $hash_table_hsh->[1]);
	$hash_table_hsh->[2] = undef unless ($hash_table_hsh->[2]);
	$hash_table_hsh->[3] = undef unless ($hash_table_hsh->[3]);
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

	&recursive_hash($SYMBOL, $syn_words_ary, $new_hsh, $synonym);

}; # sub recursive_hash


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

	my $file           = shift;
	my $stop_words_hsh = shift;
	
	open (WORDS, $file) 
		or die "Can't open stopwords file\n";
	
	while (<WORDS>) {
		
		chomp;
		
		next if (/^\s+/g); # skip blank lines
		
		my $word_comma    = "$_,";
		my $word_fullstop = "$_.";

		$stop_words_hsh->{$_}             = undef if (!exists $stop_words_hsh->{$_});
		$stop_words_hsh->{ucfirst$_}      = undef if (!exists $stop_words_hsh->{ucfirst$_} and length$_ > 1);
		$stop_words_hsh->{$word_comma}    = undef if (!exists $stop_words_hsh->{$word_comma});
		$stop_words_hsh->{$word_fullstop} = undef if (!exists $stop_words_hsh->{$word_fullstop});

	}; # while <WORDS>

}; # sub stop_words_reader


#********************************************************************************
# text_analyzer()
#
# Arguments: file to analyze
#			 synonyms hash table
#			 stopwords hash
#			 
# Returns:   array with symbols appended
#			 it finds genes in a text file
#			 it calls recursive_search()
#
#

sub text_analyzer($$) {
	
	my $file           = shift;
	my $hash           = shift;
	my $stop_words_hsh = shift;
	my $outfile        = shift;
	my @tagged_lines   = ();

	open (TEXTFILE, $file)
		or die "Can't open $file : $!\n";

	while (my $line = <TEXTFILE>) {

		chomp $line;
		my @allwords       = split /\s/, $line; 
		my @words          = grep { !exists $stop_words_hsh->{$_} } @allwords; # Removes stopwords for processing
		my @words_copy     = @words;
		my $positive_lines = "";

		for (my $i = 0; $i < @words;) {
			
			my $complete_gene = "";
			&recursive_search($words[$i], \$complete_gene, $hash, \$line, \$positive_lines, @words_copy);
						
			# skip complete gene from iteration (to avoid internal matches)
			
			if ($complete_gene) {

				my @array_count = split /\s/, $complete_gene;
				my $count       = @array_count;
				
				$i += $count; 
				splice @words_copy, 0, $count;
				
			} else {

				$i++;
				splice @words_copy, 0, 1;

			} # if

			$complete_gene = "";
			
		} # foreach word

		$positive_lines  =~ s/^\s//;
		$positive_lines .= $line if ($line);
		push @tagged_lines, $positive_lines;
		
	} # while <TEXTFILE>

	close (TEXTFILE);

	return (@tagged_lines);

}; # sub text_analyzer


#********************************************************************************
# recursive_search()
#
# Arguments: first word to analyze
#			 ref to scalar with gene name (it becomes longer as the function calls itself)
#			 hash with synonyms
#			 line that it is processing
#			 a copy of all the words in the sentence (from the first word onwards)
# 
#			 
# Returns:   nothing 
#			 it tags gene names and appends HUGO approved symbol to the sentence
#
#

sub recursive_search {
	
	my $word           = shift;
	my $complete_gene  = shift;
	my $hash           = shift;
	my $line           = shift;
	my $positive_lines = shift;
	my @words_copy     = @_;

	return unless ($word); # ends function if it runs out of words

	my $possible_gene = uc($word);
	$possible_gene    =~ s/[^A-Z0-9\s]//gi;

	if (exists $hash->{$possible_gene}) {

		if ($hash->{$possible_gene}->[0] == 0) { # if gene name has one word or else

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} 
			
			&write_line($possible_gene, $line, $complete_gene, $positive_lines, $hash);

			return;
	
		} elsif ($hash->{$possible_gene}->[0] == 1) { # if gene name has one word or else

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} 

			splice @words_copy, 0, 1;
			my $new_hash = $hash->{$possible_gene}->[1];
			&recursive_search($words_copy[0], $complete_gene, $new_hash, $line, $positive_lines, @words_copy);

		} elsif ($hash->{$possible_gene}->[0] == 2) {

			if ($$complete_gene) {

				$$complete_gene .= " " . $word;

			} else {

				$$complete_gene = $word;

			} # if gene name has one word or else

			my $new_hash = $hash->{$possible_gene}->[1];
			splice @words_copy, 0, 1;
			my $next_word = "";
			
			if ($words_copy[0]) {
				
				$next_word = uc($words_copy[0]);
				$next_word =~ s/[^A-Z0-9\s]//gi;
			
			}

			if ($words_copy[0] and exists $new_hash->{$next_word}) {

				&recursive_search($words_copy[0], $complete_gene, $new_hash, $line, $positive_lines, @words_copy);

			} else {

				&write_line($possible_gene, $line, $complete_gene, $positive_lines, $hash);

				return;
				
			} # if next word exists in hash table

		}; # if index = (0 or 1 or 2)

	return;

	} # if primary key exists

}; # sub recursive_search


#********************************************************************************
# write_line
#
# Arguments: ref to possible gene
#			 line being analyzed
#			 complete gene name found
#		     ref to var with line that will be created
#  			 ref to hash with gene symbol
# 			 
#			 
# Returns: nothing
#		   it creates output line
#
#

sub write_line {
	
	my $possible_gene  = shift;
	my $line           = shift;
	my $complete_gene  = shift;
	my $positive_lines = shift;
	my $hash           = shift;
	my $quoted_gene    = quotemeta($$complete_gene);

	$$line =~ s/^(.*?)$quoted_gene//;
				
	if ($1) {
				
		$$positive_lines .= $1 . "#$$complete_gene#&&$hash->{$possible_gene}->[2]&&";
				
	} else {
				
		$$positive_lines .= " " . "#$$complete_gene#&&$hash->{$possible_gene}->[2]&&";
				
	} # if

	return; 

} # sub write_line


#********************************************************************************
# tagged_lines_filter
#
# Arguments: array ref with tagged lines
#			 file output name 
#			 
# Returns: nothing
#		   it prints in output file those lines with one or more tags
#
#

sub tagged_lines_filter {
	
	my $tagged_lines = shift;
	my $outfile      = shift;
	
	open (OUT, "> $outfile")
		or die "Can't open $outfile : $!\n";

	foreach my $line (@$tagged_lines) {

		print OUT $line, "\n" if $line =~ m/#.+#/g ;

	} # foreach line

	close (OUT);

	return;

} # sub tagged_lines_filter