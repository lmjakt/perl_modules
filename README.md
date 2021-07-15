# Perl modules

My collection of perl modules. There are almost certainly published
perl modules that provide the same functionality as the modules
included here. However, I am somewhat too lazy to go looking for
modules that are not so difficult to write; often it takes longer to
work out how to use them than it does to write something that fits
one's own usage case better.

These modules have only been tested for my own use-cases.

## gb_download.pm
Provides a single function:

```perl
download_records( search_term, db, ret_type, retmode );
```

Downloads records associated with the `search_term` and `db`.
The `search_term` is the only mandatory argument; the others
default to, `nucleotide`, `gb` and `text`.

I have used this function to download genbank entries specified by accession
and it seems to work for that. I have not tested it for any other
purpose.

Relies on the `LWP` package.

## gb_parser.pm
Provides functions to parse and translate sequences in the genbank format.
The current version does not parse the header field correctly and doesn't
handle repeats of major fields that are often found (eg. REFERENCE). But it
does allow the extraction of all features. Handles gb files containing
multiple sequences (delimited by `//`).

Exports:

### `parse_from_file(file, is_text)`  
Reads and parses records in `file`; `is_text` is an optional parameter; if
specified and `TRUE` then `parse_from_file` will parse the `file` string
itself. Returns an array of hashes with the following structure:

- description
  - LOCUS
  - DEFINITION
  - etc...
- features
  - CDS
    - CDS_1 seq and information
	- CDS_2 seq and information
  - Gene
	- Gene 1 and information
	- Gene 2 and information
  - etc...
- sequence

Note that each field in the description is actually a hash with a key of `1`
if is the base field (eg. SOURCE) and a key of any optional sub-field
(eg. SOURCE -> ORGANISM). This doesn't support repeats of either the primary
or secondary fields.

The value of the `features` entry is an array where each feature is a hash
where the keys are either the qualifiers (eg. `/gene="AXL2"`) or the term
`range` used for the range information. The values are the unquoted text
following the `=`.

### `extract_feature_seq(feat_r, seq_r, type)`

Extracts the sequences associated with a given feature type (eg. CDS) and a
given genbank entry (i.e. an element of the array returned by
`parse_from_file`). Both the genbank entry and the sequence should be passes
as references. Returns an array of elements where each element is a hash with
two keys: `seq` and `meta`. The `meta` element is a further hash that contains
the information associated with the feature with each key being one of the
genbank feature qualifier tags. If the `codon_start` qualifier is present the
function will trim the sequence so that the first nucleotide of the sequence
corresponds to the first nucleotide of the first codon in the feature.

### `translate_peptide(seq, frame, code)`

Translates `seq` starting from `frame` using the genetic code specified by
`code`. Genetic codes are identified using the system given
at: <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>. For example,
the standard code is "1", The vertebrate mitocondrial code is "2" and so
on. Note that although numbers are used, these should be considered as strings
as they do not form a continuous sequence.

### `genetic_code(code)`

Returns a hash containing the specified translation code and a
description. The code is returned as a hash with two keys: `description`
holding a descriptor and `code` holding a hash of the code itself. If `code`
is not specified or unknown, the function will return the standard code. If
`code` is negative a hash of all functions and their descriptions will be
returned.

## `test_gb_parser.pl`
A perl script testing some of the functions in the above modules.

## genetic_codes.txt
A set of translation codes obtained from NCBI.

## `colIIA.gb` and `gb_seq_data.gb`
Two genbank formatted files. The latter one contains multiple entries.
