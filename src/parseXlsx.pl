#!/usr/bin/env perl

use strict;
use Spreadsheet::XLSX;

my $xlsx = Spreadsheet::XLSX -> new($ARGV[0]);
my $n    = 0;

foreach my $sheet (@{$xlsx -> {Worksheet}}) 
{
    $n++;
    $sheet -> {MaxRow} ||= $sheet -> {MinRow};
    foreach my $row (1 .. $sheet -> {MaxRow}) 
    {
        my $sample = $sheet->{Cells}[$row][2]->{Val};
        my $pos    = $sheet->{Cells}[$row][1]->{Val};
        if ($sheet->{Cells}[$row][8])
        {
            print "$sample,",$sheet->{Cells}[$row][8]->{Val},"\n";
        }
        else
        {
            print "$sample,${sample}_${n}_$pos\n";
        }
    }
}

