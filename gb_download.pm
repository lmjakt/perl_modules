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
my %ncbi;

## searches, and downloads records
sub download_records {
    my($search_term, $db, $ret_type, $retmode, $retmax, $api_key) = @_;
    ## this sets arguments that have not been set:
    ## // is similar to ||, but can be used to test if something is defined.
    $db //= "nucleotide";
    $ret_type //= "gb";
    $retmode //= "text";
    $retmax //= 100;
    if(!defined($api_key)){
	$api_key = look_for_api_key();
    }
    ## there are probably other modifications which should be done, but:
    ## remove spaces:
    $search_term =~ s/\s/%20/g;
    my $search_url = $base_url.$search.$db."&term=$search_term&usehistory=y";
    $search_url .= "&api_key=$api_key" if(defined($api_key));
    my $web_env_string = get($search_url);
    my( $web, $key, $count ) = extract_key( $web_env_string );
    my $data = "";
    for( my $retstart=0; $retstart < $count; $retstart += $retmax ){
	my $fetch_url = $base_url.$fetch.$db."&query_key=$key&WebEnv=$web"."&rettype=$ret_type&retmode=$retmode".
	    "&retstart=$retstart&retmax=$retmax";
	$fetch_url .= "&api_key=$api_key" if(defined($api_key));
	$data .= get($fetch_url);
    }
    return($data);
}

sub extract_key {
    my $str = shift @_;
    my $web = $1 if ($str =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($str =~ /<QueryKey>(\d+)<\/QueryKey>/);
    my $count = $1 if ($str =~ /<Count>(\d+)<\/Count>/);
    return($web, $key, $count);
}
 
## returns an undefined value if it finds nothing.   
sub look_for_api_key {
    my $key;
    return($ENV{NCBI_API_KEY}) if(defined($ENV{NCBI_API_KEY}));
    read_config() if(!%ncbi);
    return($ncbi{api_key}) if(defined($ncbi{api_key}));
    return($key);
}

sub read_config {
    my $home = $ENV{HOME};
    my @conf_files = ("$home/.ncbi.conf", "$home/.ncbi/ncbi.conf");
    for my $file(@conf_files){
	if( -r $file ){
	    open(my $in, "<", $file) || next;
	    while(<$in>){
		chomp;
		my @tmp = split /=/, $_;
		if(@tmp > 1){
		    $ncbi{$tmp[0]} = $tmp[1];
		}
	    }
	    last;
	}
    }
}

1;
