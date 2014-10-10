#!/usr/bin/perl
#
# gene_finder.pl -
#	This script reads a file that contains a list of gene symbols
#	and its corresponding synonyms and it creates a synonyms hash table.
#	Then, it looks for matches in a text file and saves the lines
#	containing them as matches_originalfile.txt
#
#	Arguments:	Synonyms tabular file
#				Stopwords file
#				Text file(s) to analyze
#
#
#********************************************************************************


use warnings;
use strict;
use Data::Dumper;
use LWP::Simple;
use utf8;

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

table_reader($synonyms);
my %greek_dict = init_greek();

print STDERR "# READING STOPWORDS LIST...\n\n";

stop_words_reader($stop_words_file, \%stop_words);

print STDERR "# LOOKING FOR GENES IN TEXT...\n";

foreach my $file (@problem_files) {
	print STDERR "\t# Analyzing $file...\n";

	my @tagged_lines = text_analyzer($file, \%hash_table, \%stop_words, \%greek_dict);

	$file =~ s/.+\///;
	my $matches_out  = "matches_$file"; # Creating output file name
	tagged_lines_filter(\@tagged_lines, $matches_out);

	print STDERR "\t# Matches saved as: $matches_out\n\n";
}

print STDERR "\n## PROGRAM FINISHED ##\n";


#********************************************************************************
# FUNCTIONS
#********************************************************************************

#********************************************************************************
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

		next if (length $Ap_SYMBOL < 3);

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
				$match    =~ s/^\s//;
    			push @quoted_fields, $match;
   			 }; # while match
		};

		@all_synonyms = (@normal_fields, @quoted_fields);

		@all_synonyms = grep { length($_) > 2} @all_synonyms;
			# remove short synonyms

		@all_synonyms = grep { $_ =~ m/[a-z]/gi } @all_synonyms;
			# remove numeric synonyms

		remove_synonyms(\@all_synonyms);
		hash_creator($Ap_SYMBOL, \@all_synonyms);
	} # while <SYN>

	close (SYN);

} # sub table_reader


#********************************************************************************
sub init_greek {

	my %greek_dict = (
					α => 'alpha',
					β => 'beta',
					γ => 'gamma',
					δ => 'delta',
					ε => 'epsilon',
					ζ => 'zeta',
					η => 'eta',
					θ => 'theta',
					ι => 'iota',
					κ => 'kappa',
					λ => 'lambda',
					μ => 'mu',
					ν => 'nu',
					ξ => 'xi',
					π => 'pi',
					ρ => 'rho',
					σ => 'sigma',
					τ => 'tau',
					τ => 'upsilon',
					φ => 'phi',
					χ => 'chi',
					ψ => 'psi',
					ω => 'omega'
					);

	return (%greek_dict);

} # sub greek_dict


#********************************************************************************
sub remove_synonyms {

	my $synonyms_array = shift;
	my @not_valid      = ( 'SARA',
						   'ROS',
						   'P200',
						   'P1.10',
						   'DAG',
						   'gamma',
						   'beta',
						   'delta',
						   'anova',
						   'GMP',
						   'RNA binding protein',
						   'RNA-binding protein',
						   'neuronal migration',
						   'yes',
						   'minor',
						   'p110',
						   'proc',
						   'mass',
						   'red',
						   'BETA-3',
						   'Beta3',
						   'beta2',
						   'mice',
						   'cell',
						   'rho',
						   'HCC',
						   'LPS',
						   'chip');

	foreach my $stopword (@not_valid) {
		@$synonyms_array = grep { $_ !~ m/^$stopword$/gi} @$synonyms_array;
	}

	return;

}

#********************************************************************************
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
			recursive_hash($SYMBOL, \@syn_words, $hash_table{$first_word}, $synonym);
	    } # if

	} # foreach

} # sub hash_creator


#********************************************************************************
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

	recursive_hash($SYMBOL, $syn_words_ary, $new_hsh, $synonym);

}; # sub recursive_hash


#********************************************************************************
sub stop_words_reader {

	my $file           = shift;
	my $stop_words_hsh = shift;

	open (WORDS, $file)
		or die "Can't open stopwords file\n";

	while (<WORDS>) {
		chomp;
		next if (/^\s+/g); # skip blank lines

		my @forms = ();
		push @forms, "$_";
		push @forms, "$_,";
		push @forms, "$_.";
		push @forms, "$_;";
		push @forms, "($_";
		push @forms, "$_)";

		foreach my $form (@forms) {
			$stop_words_hsh->{$form} = undef;
			$stop_words_hsh->{ucfirst$form} = undef if (length$_ > 1 and substr($form, 0, 1) ne "(");
		}

	}; # while <WORDS>

}; # sub stop_words_reader


#********************************************************************************
sub text_analyzer {

	my $file           = shift;
	my $hash           = shift;
	my $stop_words_hsh = shift;
	my $greek_dict     = shift;
	my @tagged_lines   = ();

	open TEXTFILE, '<:encoding(UTF-8)', $file
		or die "Can't open $file : $!\n";

	while (my $line = <TEXTFILE>) {

		chomp $line;
		$line =~ s/(\p{Greek})/$greek_dict->{lc$1}/ge; # Greek letters to ASCII
		my @allwords       = split /\s/, $line;
		my @words          = grep { !exists $stop_words_hsh->{$_} } @allwords; # Removes stopwords for processing
		my @words_copy     = @words;
		my $positive_lines = "";

		for (my $i = 0; $i < @words;) {
			my $complete_gene = "";
			recursive_search($words[$i], \$complete_gene, $hash, \$line, \$positive_lines, @words_copy);

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

			write_line($possible_gene, $line, $complete_gene, $positive_lines, $hash);

			return;

		} elsif ($hash->{$possible_gene}->[0] == 1) { # if gene name has one word or else

			if ($$complete_gene) {
				$$complete_gene .= " " . $word;
			} else {
				$$complete_gene = $word;
			}

			splice @words_copy, 0, 1;
			my $new_hash = $hash->{$possible_gene}->[1];
			recursive_search($words_copy[0], $complete_gene, $new_hash, $line, $positive_lines, @words_copy);
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
				recursive_search($words_copy[0], $complete_gene, $new_hash, $line, $positive_lines, @words_copy);
			} else {
				write_line($possible_gene, $line, $complete_gene, $positive_lines, $hash);
				return;
			} # if next word exists in hash table

		}; # if index = (0 or 1 or 2)

	return;

	} # if primary key exists

}; # sub recursive_search


#********************************************************************************
sub write_line {

	my $possible_gene  = shift;
	my $line           = shift;
	my $complete_gene  = shift;
	my $positive_lines = shift;
	my $hash           = shift;

	$$complete_gene =~ s/[\.\),;]+$//;	# Word boundaries have to be removed
	$$complete_gene =~ s/^[\.\(,]//;	# so they can be added in the regex

	my $quoted_gene = $$complete_gene;
	$quoted_gene    = quotemeta($quoted_gene);
	my $bound = '(?:(?<![+\w/\-\p{Greek}])(?=[+\w/\-\p{Greek}])|(?<=[+\w/\-\p{Greek}])(?![+\w/\-\p{Greek}]))';

	my $flag = $$line =~ s/^(.*?)$bound$quoted_gene$bound//;
						# Word boundaries are necessary so the
						# tagging will be done in genes and not
						# in words that contain gene names
						# eg. barMAPK, ERKfoo...
	if ($flag) {

		if ($1) {
			$$positive_lines .= $1 . "#$$complete_gene#&&$hash->{$possible_gene}->[2]&&";
		} else {
			$$positive_lines .= " " . "#$$complete_gene#&&$hash->{$possible_gene}->[2]&&";
		} # if words before gene

	} # if match

	return;

} # sub write_line


#********************************************************************************
sub tagged_lines_filter {

	my $tagged_lines = shift;
	my $outfile      = shift;

	open (OUT, '>:utf8', $outfile)
		or die "Can't open $outfile : $!\n";

	foreach my $line (@$tagged_lines) {
		print OUT $line, "\n" if $line =~ m/#.+#.+#.+#/g ;
	} # foreach line

	close (OUT);
	return;

} # sub tagged_lines_filter
