package gb_parser;
require(Exporter);
@ISA = qw(Exporter);
@EXPORT = qw( parse_from_file extract_feature_seq translate_peptide genetic_code );

use warnings;
use strict;

## BUGS
## Currently only returns a single REFERENCE section
## This is problematic since the REFERENCE field may
## be repeated several times in the header part.
## Modifying this requires returning a
## hash of arrays of hashes rather than a simple
## hash of hashes.


my %functions = (
    "join" => \&parse_operations,
    "order" => \&parse_operations,
    "complement" => \&complement
    );

sub description_fields {
    return( qw( LOCUS DESCRIPTION ACCESSION VERSION KEYWORDS SOURCE REFERENCE ) );
}

## argument can either be a string containing a file name
## or a string containing the actual data. If the latter, then
## $is_text should be set to tru.
sub parse_from_file {
    my ($file, $is_text) = @_;
    my @lines = ();
    my @entries = ();
    $is_text = 0 if !defined($is_text);
    my $in;
    if(!$is_text){
	open($in, "<", $file) || die "unable to open $file $!\n";
    }else{
	open($in, "<", \$file) || die "unable to open in memory file $!\n";
    }
    while( @lines = read_entry($in) ){
	my %description = ();
	my %features = ();
	my $sequence = "";
	my $i = 0;
	while($i < @lines){
	    while( ($i = read_description(\@lines, $i, \%description)) ){
		last if( defined($description{"FEATURES"}) || $i < 0 );
	    }
	    $i = read_features(\@lines, $i, \%features);
	    $sequence = read_sequence(\@lines, $i);
	    last;
	}
	push @entries, {"description"=>{%description}, 
			    "features"=>{%features},
			    "sequence"=>$sequence };
    }
    return(@entries);
}

    
sub read_entry {
    my($in) = @_;
    my @lines = ();
    while(<$in>){
	chomp;
	last if( $_ =~ m#^//# );
	push @lines, $_ if $_ =~ /\S+/;
    }
    return(@lines);
}

sub read_description {
    my($lines, $i, $des_ref) = @_;
    my $key = "";
    my $value_ref = 0;
    if($lines->[$i] =~ /^(\S+)\s+(.+)/){
	$key = $1;
	$des_ref->{$key}{1} = $2;
	return( $i + 1 ) if $key eq "FEATURES";
	$value_ref = \$des_ref->{$key}{1};
    }else{
	return(-1);
    }
    $i++;
    while($i < @{$lines} && $lines->[$i] =~ /^\s+(.+)/){
	if($lines->[$i] =~ /^\s{12}(\S+.+)/){
		${$value_ref} .= " ".$1;
	}
	if($lines->[$i] =~ /^\s{2}(\S+)\s+(.+)/){
	    $des_ref->{$key}{$1} = $2;
	    $value_ref = \$des_ref->{$key}{$1};
	}
	$i++;
    }
    return($i);
}

sub unquote {
    my $ref = shift @_;
    return if(!$ref);
    ${$ref} =~ s/^\s*"//;
    ${$ref} =~ s/"\s*$//;
}

sub read_features {
    my($lines, $i, $feat_ref) = @_;
    my $value_ref = 0;
    my $key = "";
    my $sep = "";
    my $j;
    for($j=$i; $j < @{$lines}; $j++){
	last if( $lines->[$j] =~ /^ORIGIN/ );
	if( $lines->[$j] =~ /^\s{5}(\S+)\s+(.+)/ ){
	    $key = $1;
	    my $value = $2;
	    unquote( $value_ref );
	    push @{$feat_ref->{$key}}, { range => $value };
	    $value_ref = \$feat_ref->{$key}[-1]{range};
	    $sep = "";
	    next;
	}
	if( $lines->[$j] =~ m#^\s{21}/([^=]+)=(.+)# ){
	    my $key2 = $1;
	    my $value = $2;
	    unquote( $value_ref );
	    $feat_ref->{$key}[-1]{$key2} = $value;
	    $value_ref = \$feat_ref->{$key}[-1]{$key2};
	    $sep = $key2 eq "translation" ? "" : " ";
	    next;
	}
	if( $lines->[$j] =~ /^\s{21}(\S+.+)/ ){
	    ${$value_ref} .= ($sep.$1);
	}
    }
    unquote($value_ref);
    return($j);
}

sub extract_feature_seq {
    my($feat_r, $seq_r, $type) = @_;
    return(0) if !defined($feat_r->{$type});
    my @seq_data = ();
    for my $hash_r( @{$feat_r->{$type}} ){
	my $range_term = $hash_r->{range};
	$range_term =~ s/\s//g;
	my $sequence = parse_operations($range_term, $$seq_r);
	if(defined($hash_r->{codon_start})){
	    $sequence = substr($sequence, $hash_r->{codon_start}-1);
	}
	push @seq_data, {seq => $sequence, meta => {%{$hash_r}} };
    }
    return @seq_data;
}

sub read_sequence {
    my($lines, $i) = @_;
    my $sequence = "";
    $i++ if($lines->[$i] =~ /ORIGIN/);
    while($i < @{$lines} && $lines->[$i] !~ m#^//#){
	chomp($lines->[$i]);
	$lines->[$i] =~ s/[\s0-9]//g;
	$sequence .= $lines->[$i];
	$i++;
    }
    return($sequence);
}
## should only take a single range,
sub extract_range {
    my($range, $string) = @_;
    if($range =~ /^(\d+)&/){
	return(substr( $string, $1, 1 ));
    }
    if($range =~ /^<?(\d+)\.\.>?(\d+)>?$/){
	return(substr( $string, $1-1, 1 + $2-$1 ));
    }
    die "Unexpected range term: $range\n";    
}

sub complement {
    my $arg_string = shift @_;
    my $seq = shift @_;
    $arg_string = extract_br_contents($arg_string);
    my $string = parse_operations($arg_string, $seq);
    ## reverse complement and then return
    $string = reverse($string);
    $string =~ tr/ACGTUMRWSYKVHDBN/TGCAAKYWSRMBDHVN/;
    $string =~ tr/acgtumrwsykvhdbn/tgcaakywsrmbdhvn/;
    return($string);
}

# sub join {
#     my $arg_string = shift @_;
#     my $seq = shift @_;
#     ## extract the content of the operation
#     return( parse_operations($arg_string, $seq) );
# }
    
sub parse_operations {
    my $expression = shift @_;
    my $seq = shift @_;
    my $string = "";
    while(length($expression)){
	my($first, $remainder) = extract_argument( $expression );
	if($first =~ /([^(]+)\((.+)\)$/){
	    $string .= $functions{$1}( $2, $seq );
	}else{
	    $string .= extract_range( $first, $seq );
	}
	$expression = $remainder;
    }
    return($string);
}

sub extract_argument {
    my $arg_string = shift @_;
    my $comma_pos = index( $arg_string, "," );
    my $br_pos = index($arg_string, "(");
    return( $arg_string, "" ) if($comma_pos == -1 && $br_pos == -1);
    goto COMMA if( $br_pos == -1 );
    goto BRACKET if( $comma_pos == -1 );
    goto COMMA if( $comma_pos < $br_pos );
    goto BRACKET;

  COMMA:
    return( substr($arg_string, 0, $comma_pos), 
	    substr($arg_string, $comma_pos + 1 ));
  BRACKET:
    my $i;
    my $l_count = 0;
    my $r_count = 0;
    for($i=$br_pos; $i < length($arg_string); ++$i){
	$l_count++ if substr($arg_string, $i, 1) eq "(";
	$r_count++ if substr($arg_string, $i, 1) eq ")";
	last if $l_count == $r_count;
    }
    my $extra = substr($arg_string, $i + 1, 1) eq "," ? 1 : 0;
    return( substr($arg_string, 0, $i+1),
	    substr($arg_string, $i + 1 + $extra));
}

sub extract_br_contents {
    my $op = shift @_;
    $op =~ s/^[^(]+\(//;
    $op =~ s/\($//;
    return($op);
}

## frame is given in 1-based coordinates since this is what
## is present in the gb file. Hence the decrement at the beginning
## of the function.
sub translate_peptide {
    my($seq, $frame, $code) = @_;
    $seq = uc($seq);
    $frame--;
    print "frame is $frame\n";
    my $aa = "";
    my %gene_code = genetic_code($code);
    my %gcode = %{$gene_code{code}};

    for(my $i=$frame; $i < length($seq)-2; $i += 3){
	if(defined($gcode{substr($seq, $i, 3)})){
	    $aa .= $gcode{substr($seq, $i, 3)};
	}else{
	    $aa .= "X";
	}
    }
    print "\n";
    return($aa);
}

## give the genetic code id (a number treated as a character)
## as given on:
## https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
## if no code id given, return the standard code ("1")
## if a negative numer is given return all codes.
sub genetic_code {
    my $code = shift @_;
    if(!defined($code)){
	$code = "1"
    }
    my %gcode;
    $gcode{"1"}{description} = "The Standard Code";
    %{$gcode{"1"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"2"}{description} = "The Vertebrate Mitochondrial Code (transl_table=2)";
    %{$gcode{"2"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"M",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"*",  "AGG"=>"*",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"3"}{description} = "The Yeast Mitochondrial Code (transl_table=3)";
    %{$gcode{"3"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"T",  "CTC"=>"T",
			     "CTA"=>"T",  "CTG"=>"T",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"M",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"4"}{description} = "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)";
    %{$gcode{"4"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"5"}{description} = "The Invertebrate Mitochondrial Code (transl_table=5)";
    %{$gcode{"5"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"M",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"S",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"6"}{description} = "The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)";
    %{$gcode{"6"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Q",  "TAG"=>"Q",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"9"}{description} = "The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)";
    %{$gcode{"9"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			     "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			     "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			     "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			     "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			     "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			     "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			     "AAA"=>"N",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"S",
			     "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			     "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			     "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"10"}{description} = "The Euplotid Nuclear Code (transl_table=10)";
    %{$gcode{"10"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"C",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"11"}{description} = "The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)";
    %{$gcode{"11"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"12"}{description} = "The Alternative Yeast Nuclear Code (transl_table=12)";
    %{$gcode{"12"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"S",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"13"}{description} = "The Ascidian Mitochondrial Code (transl_table=13)";
    %{$gcode{"13"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"M",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"G",  "AGG"=>"G",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"14"}{description} = "The Alternative Flatworm Mitochondrial Code (transl_table=14)";
    %{$gcode{"14"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Y",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"N",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"S",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"16"}{description} = "Chlorophycean Mitochondrial Code (transl_table=16)";
    %{$gcode{"16"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"L",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"21"}{description} = "Trematode Mitochondrial Code (transl_table=21)";
    %{$gcode{"21"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"M",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"N",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"S",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"22"}{description} = "Scenedesmus obliquus Mitochondrial Code (transl_table=22)";
    %{$gcode{"22"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"*",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"L",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"23"}{description} = "Thraustochytrium Mitochondrial Code (transl_table=23)";
    %{$gcode{"23"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"*",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"24"}{description} = "Rhabdopleuridae Mitochondrial Code (transl_table=24)";
    %{$gcode{"24"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"K",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"25"}{description} = "Candidate Division SR1 and Gracilibacteria Code (transl_table=25)";
    %{$gcode{"25"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"G",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"26"}{description} = "Pachysolen tannophilus Nuclear Code (transl_table=26)";
    %{$gcode{"26"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"*",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"A",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"27"}{description} = "Karyorelict Nuclear Code (transl_table=27)";
    %{$gcode{"27"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Q",  "TAG"=>"Q",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"28"}{description} = "Condylostoma Nuclear Code (transl_table=28)";
    %{$gcode{"28"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Q",  "TAG"=>"Q",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"29"}{description} = "Mesodinium Nuclear Code (transl_table=29)";
    %{$gcode{"29"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Y",  "TAG"=>"Y",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"30"}{description} = "Peritrich Nuclear Code (transl_table=30)";
    %{$gcode{"30"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"E",  "TAG"=>"E",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"*",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"31"}{description} = "Blastocrithidia Nuclear Code (transl_table=31)";
    %{$gcode{"31"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"E",  "TAG"=>"E",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"R",  "AGG"=>"R",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    $gcode{"33"}{description} = "Cephalodiscidae Mitochondrial UAA-Tyr Code (transl_table=33)";
    %{$gcode{"33"}{code}} = ( "TTT"=>"F",  "TTC"=>"F",  "TTA"=>"L",  "TTG"=>"L",  "TCT"=>"S",  "TCC"=>"S",
			      "TCA"=>"S",  "TCG"=>"S",  "TAT"=>"Y",  "TAC"=>"Y",  "TAA"=>"Y",  "TAG"=>"*",
			      "TGT"=>"C",  "TGC"=>"C",  "TGA"=>"W",  "TGG"=>"W",  "CTT"=>"L",  "CTC"=>"L",
			      "CTA"=>"L",  "CTG"=>"L",  "CCT"=>"P",  "CCC"=>"P",  "CCA"=>"P",  "CCG"=>"P",
			      "CAT"=>"H",  "CAC"=>"H",  "CAA"=>"Q",  "CAG"=>"Q",  "CGT"=>"R",  "CGC"=>"R",
			      "CGA"=>"R",  "CGG"=>"R",  "ATT"=>"I",  "ATC"=>"I",  "ATA"=>"I",  "ATG"=>"M",
			      "ACT"=>"T",  "ACC"=>"T",  "ACA"=>"T",  "ACG"=>"T",  "AAT"=>"N",  "AAC"=>"N",
			      "AAA"=>"K",  "AAG"=>"K",  "AGT"=>"S",  "AGC"=>"S",  "AGA"=>"S",  "AGG"=>"K",
			      "GTT"=>"V",  "GTC"=>"V",  "GTA"=>"V",  "GTG"=>"V",  "GCT"=>"A",  "GCC"=>"A",
			      "GCA"=>"A",  "GCG"=>"A",  "GAT"=>"D",  "GAC"=>"D",  "GAA"=>"E",  "GAG"=>"E",
			      "GGT"=>"G",  "GGC"=>"G",  "GGA"=>"G",  "GGG"=>"G" );
    
    if($code < 0){
	return(%gcode);
    }
    if(!defined($gcode{$code})){
	print STDERR "No code defined for $code\nWill use standard code instead\n";
	$code = "1";
    }
    return( %{$gcode{$code}} );
}

1;
   
