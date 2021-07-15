#!/usr/bin/perl -w
BEGIN {
    push @INC, "./";
}
use warnings;
use strict;
use gb_parser;
use gb_download;

my $file = shift @ARGV;

my @seq_data = parse_from_file( $file );

print STDERR "obtained info from ", scalar(@seq_data), " sequences\n";

for my $hash_ref(@seq_data){
    print $hash_ref->{description}{LOCUS}{1}, "\n";
#   ## extract all CDSs and translate with the standard table.
    my @cds_data = extract_feature_seq( $hash_ref->{features}, \$hash_ref->{sequence}, "CDS" );
    next if(!$cds_data[0]);
    for my $seq_dr(@cds_data){
	print "range:\n$seq_dr->{meta}{range}\n";
	print "sequence:\n", $seq_dr->{seq}, "\n";
	print "translation:\n$seq_dr->{meta}{translation}\n";
	print translate_peptide( $seq_dr->{seq}, 1 ), "\n";
    }
}

## make use of gb_download to download these sequences

my $seq_id = "XM_011700910,KM273030,AY310175,XM_017577217,XM_008911095,XM_004670360,LC102263.1,KT726962.1,HM237355.1";

my @ids = split /,/, $seq_id;
my $data = download_records( $seq_id );

for my $hash_ref(@seq_data){
    print "\n\n", $hash_ref->{description}{LOCUS}{1}, "\n";
     ## extract all CDSs and translate using the specified table
    my @cds_data = extract_feature_seq( $hash_ref->{features}, \$hash_ref->{sequence}, "CDS" );
    next if(!$cds_data[0]);
    for my $seq_dr(@cds_data){
	print "range:\n$seq_dr->{meta}{range}\n";
	print "sequence:\n", $seq_dr->{seq}, "\n";
	print "translation:\n$seq_dr->{meta}{translation}\n";
	my @keys = keys %{$seq_dr->{meta}};
	print "keys: ", join(", ", @keys), "\n";
	my $trans_table = $seq_dr->{meta}{transl_table} || "1";
	print "translation table: $trans_table\n";
	print translate_peptide( $seq_dr->{seq}, 1, $trans_table), "\n";
    }
}
