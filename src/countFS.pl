#!/usr/bin/env perl

use BioSeqFormat;
use strict;

my $file = $ARGV[0];
my $fh;
open($fh,$file);
my $f  = fasta->new();
while ($f->to_object($fh))
{
    my $s = $f->get_seq();
    print $file,",";
    print $f->get_id(),",";
    ($s =~ /\*\w/)? print "FS\n" : print "OK\n";
}
close($fh);
