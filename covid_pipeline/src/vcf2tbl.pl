#!/usr/bin/env perl

use Getopt::Long;
use strict;

use vars qw($sample_name);

GetOptions(
    "n=s"   => \$sample_name
);

my $name;
while (<>)
{
    chomp;
    my @v = split /\s+/;
    
    if (!$name)
    {
        $name = $sample_name || $v[0];
    }
    
    if ($v[6] eq "PASS" || $v[6] eq '.')
    {
        my $pos = $v[1];
        my $wt  = $v[3];
        my $mt  = $v[4];
        my $sc  = $v[5];

        my $status = 'M';
        my $delta  = 0;
        if ( length($wt) > length($mt) )
        {
            $status = 'D';
            $delta  = length($wt) - length($mt);
        }
        if ( length($wt) < length($mt) )
        {
            $status = 'I';
            $delta  = length($mt) - length($wt);
        }
        
        print "$name,$pos,$status,$delta,$wt,$mt,$sc\n";
    }
}
