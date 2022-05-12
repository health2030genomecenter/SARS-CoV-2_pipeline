#!/usr/bin/env perl

use strict;
my %db;
my $fh;

# Read DIANALABS samplesheet
open($fh, $ARGV[0]) || die;
while (<$fh>)
{
    my @v = split /;/;
    $v[2] =~ s/"//g;
    $db{$v[2]} = 1;
}
close($fh);

# Read GC samplesheet
open($fh, $ARGV[1]) || die;
while (my $l = <$fh>)
{
    my @v = split /,/,$l;
    my ($t) = $v[1] =~ /^(.*)_null_*/;
    if ($v[9] =~ /DIANALABS/)
    {
        if ($db{$t})
        {
            print $l;
            warn("Keep $t\n");
        }
        else
        {
            $v[1] = "0$v[1]";
            $v[2] = "0$v[2]";
            my ($t2) = $v[1] =~ /^(.*)_null_*/;
            if ($db{$t2})
            {
                warn("Modifiy $t -> $t2\n");
                print join(",",@v);
            }
            else 
            {
                warn "WRONG SAMPLE ID $t\n";
                print $l;
            }
        }
    }
    else
    {
        print $l;
    }
}
close($fh)
#

