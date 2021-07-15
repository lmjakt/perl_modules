package gb_download;
require(Exporter);
@ISA = qw(Exporter);
@EXPORT = qw( download_records );

use warnings;
use strict;
use LWP::Simple;

## BUGS
## No error handling whatsoever implemented at the moment

my $base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $fetch = "efetch.fcgi?db=";
my $search = "esearch.fcgi?db=";
my $summary = "esummary.fcgi?db=";

## searches, and downloads records
sub download_records {
    my($search_term, $db, $ret_type, $retmode) = @_;
    ## this sets arguments that have not been set:
    ## // is similar to ||, but can be used to test if something is defined.
    $db //= "nucleotide";
    $ret_type //= "gb";
    $retmode //= "text";
    ## there are probably other modifications which should be done, but:
    ## remove spaces:
    $search_term =~ s/\s/%20/g;
    my $search_url = $base_url.$search.$db."&term=$search_term&usehistory=y";
    my $web_env_string = get($search_url);
    my( $web, $key ) = extract_key( $web_env_string );
    my $fetch_url = $base_url.$fetch.$db."&query_key=$key&WebEnv=$web"."&rettype=$ret_type&retmode=$retmode";
    my $data = get($fetch_url);
    return($data);
}

sub extract_key {
    my $str = shift @_;
    my $web = $1 if ($str =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($str =~ /<QueryKey>(\d+)<\/QueryKey>/);
    return($web, $key);
}
    

# $db="nucleotide";
# $query=$seq_id;

# $url=$base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

# $output = get($url);
# print "output:\n$output\n\n";

# ## parse the webenv and query key?
# $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
# $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

# ## summary url:
# $url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

# ## get the summary:
# $docsums = get($url);
# print $docsums;

# $url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
# $url .= "&rettype=gb&retmode=text";
# ## get more:
# $data = get($url);

# open(OUT, ">gb_seq_data.gb") || die "unable to open gb_seq_data.gb $!\n";
# print OUT "$data";
# close(OUT);
