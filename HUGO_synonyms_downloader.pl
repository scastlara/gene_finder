#!/usr/bin/perl
#
# HUGO_synonyms_downloader.pl - THis script downloads 
#

use warnings;
use strict;
use LWP::Simple;


my $outfile = "HUGO_synonyms.tbl";

my $url = 'http://www.genenames.org/cgi-bin/download?'.
          'col=gd_app_sym&'.        
          'col=gd_app_name&'.
          'col=gd_prev_sym&'.
          'col=gd_aliases&'.
          'col=gd_pub_ensembl_id&'.
          'col=gd_pub_refseq_ids&'.
          'col=md_prot_id&'.
          'col=gd_prev_name&'.
          'col=gd_name_aliases&'.
          'status=Approved&'.
          'status=Entry%20Withdrawn&'.
          'status_opt=2&'.
          'where=&'.
          'order_by=gd_hgnc_id&'.
          'format=text&'.
          'limit=&'.
          'hgnc_dbtag=on&'.
          'submit=submit';

my $page = get($url);

$page =~ s/FIG, //g;

open (OUT, ">$outfile")
	or die "Can't open $outfile : $! \n";

print OUT $page;

close (OUT);