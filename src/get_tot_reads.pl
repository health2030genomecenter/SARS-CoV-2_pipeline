#!/usr/bin/env perl

use strict;
use lib '.';
use JSON;
use Data::Dumper;

my %db;
my $fh;
if ( $ARGV[0] eq '-' )
{
    $fh = \*STDIN;
    parse($fh,\%db);
    close($fh);
}
else
{
    foreach my $f (@ARGV)
    {
        open( $fh, "<", $f ) || die "Cannot read JSON file '$f'.\n";
        parse($fh,\%db);
        close($fh);
    }
}

foreach my $k (keys %db)
{
    print "$k,",2*$db{$k},"\n"; # get the reads (2*cluster)
}


sub parse 
{
    my ($fh,$db) = @_;
    my $json_string = do { local $/; <$fh> };
    my $json = JSON::XS->new->decode($json_string);
    foreach my $a (@{$json->{ConversionResults}})
    {
        my $b = $a->{DemuxResults};
        foreach my $c (@$b)
        {
            $db{$c->{SampleId}} += $c->{NumberReads};
        }
    }
}

