#!/usr/bin/env perl

# USAGE
#  ./omicronINS22200.pl variant_attribition.csv variant_attribition.log

use Data::Dumper;
use strict;

my %db;
my $fh;

# Read variant_attribition.csv
open($fh, $ARGV[0]) || die "Missing input variant_attribution.csv file.\n";
while (<$fh>)
{
    chomp;
    my @v = split /,/;
    
    if ($v[1] =~ /Omicron/)
    {
        $db{$v[0]} = ($v[2] < 1)? 1 : 2;
    }
}
close($fh);

# Read variant_attribition.log
open($fh, $ARGV[1]) || die "Missing input variant_attribution.log file.\n";

my $buffer1 = '';
my $sample = '';

while (<$fh>)
{
    chomp;
    
    if (/^#SAMPLE: (.*)/ && $db{$1})
    {
        $sample = $1;
        next;
    }

    if (/^#\/\// && $sample)
    {
        if ($db{$sample} > 1)
        {
            print "$sample,OK\n";
        }
        else
        {
            print "$sample,$buffer1\n";
        }
        $sample = '';
        $buffer1 = '';
        next;
    }

    if ($sample && /^#POSITIONS: Omicron[^\,]+,(\d+),(\d+),obs=(\S+),exp=(\S+),.*(SCORE=(\S+)|MISSING)/)
    {
        $buffer1 .= ",$1:$4:$3" if ($5 eq 'MISSING' || $6 < 1); 
        next;
    }
    
    if ($sample && /^#SCORE: ([^,\/]+)/ && $1 ne "Omicron")
    {
        $sample = '';
        $buffer1 = '';
    }
}
close($fh);

