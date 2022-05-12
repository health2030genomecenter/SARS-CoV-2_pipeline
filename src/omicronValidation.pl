#!/usr/bin/env perl

# USAGE
#  ./omicronINS22200.pl variant_attribition.csv variant_attribition.log

use Data::Dumper;
use strict;

my %db;
my %vr;
my $fh;

# Read variant_attribition.csv
open($fh, $ARGV[0]) || die "Missing input variant_attribution.csv file.\n";
while (<$fh>)
{
    chomp;
    my @v = split /,/;
    
    if ($v[1] =~ /^Omicron/)
    {
        $db{$v[0]} = ($v[2] < 1)? 1 : 2;
        $vr{$v[0]} = $v[1]; 
    }
}
close($fh);

# Read variant_attribition.log
open($fh, $ARGV[1]) || die "Missing input variant_attribution.log file.\n";

print <<EOF;
# Table describing the status of the lineage descriminating mutations in the Spike gene. Fields are separated by commas.
#
# Column1: sample name
# Column2: status of lineage mutations. OK (all mutations are found), otherwise empty
# Column3 (and following): list of positions with anomalies. The format of the cells can be either
#  1) <POSITION>:not_covered - the position has no coverage.
#  2) <POSITION>:<EXPECTED_STATUS>:REF - the position contains the reference nucleic acid (NC_045512.2) instead of the expected mutation
#  3) <POSITION>:<EXPECTED_STATUS>:<OBSERVED_STATUS> - the position contains an unexpected nucleic acid
#
# Description of the cell format:
#  <POSITION>: position on the reference genome (NC_045512.2).
#  <EXPECTED_STATUS> : expected observation (format see below).
#  <OBSERVATION_STATUS> : actual observation (format see below)
#
# The status observation has the format (A:B), whew
#  * A=M (match) and B=nucleic acid - the position contains the specified nucleic acid.
#  * A=D (deletion) or A=I (insertion), and B=length of deletion or insertion

EOF

my $buffer = '';
my $sample = '';
my ($n,$m) = (0,0);

while (<$fh>)
{
    chomp;

    if (/^#SAMPLE: (\S+)/)
    {
        ($n,$m) = (0,0);
        $sample = ($vr{$1})? $1 : undef;
    }
    elsif (/^#\/\// && $sample)
    {
        my $s = $n/$m;
        if ($s < 1)
        {
            print "$sample,$buffer\n";
        }
        else
        {
            print "$sample,OK\n";
        }
        $sample = '';
        $buffer = '';
    }
    elsif ($sample && /^#POSITIONS: $vr{$sample},(\d+),(\d+),obs=([^,]+),exp=([^,]+),.*\[S\].*SCORE=(\S+)/)
    {
        $m++;
        $n++;
    }
    elsif ($sample && /^#POSITIONS: $vr{$sample},(\d+),(\d+),obs=([^,]+),exp=([^,]+).*\[S\].*(MISSING_MUTATION)/)
    {
        $m++;
        $buffer .= ",$1:$4:REF"; 
    }
    elsif ($sample && /^#POSITIONS: $vr{$sample},(\d+),(\d+),obs=([^,]+),exp=([^,]+).*\[S\].*NO_COVERAGE/)
    {
        $m++;
        $buffer .= ",$1:not_covered"; 
    }
    elsif ($sample && /^#SCORE: ([^,]+)/ && $1 ne $vr{$sample})
    {
        ($m,$n) = (0,0);
        $sample = '';
        $buffer = '';
    }
}
close($fh);

