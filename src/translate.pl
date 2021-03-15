#!/usr/bin/env perl
#
use lib "perl5lib";

use vars qw($opt_f $opt_h $opt_g $opt_x);
use Getopt::Std;
use BioSeqFormat;
use strict;

$opt_f = -1000;
getopts('f:gxh');

if ($opt_h) {

print <<EOF;

translate.pl 
    Translate DNA sequences into peptides.

SYNOPSIS:
    translate.pl [-f <s>] [-g] [-h]  DNA_FASTA

OPTIONS:
  -f <s>:   Specify frame -1,-2,-3,1,2,3, or '+' (the 3 frames in the positive 
            strand) or '-' (the 3 frames in the negative strand). All the 6 
            frames are explored if the -f option is not specified.
  -g:       Guess best frame and return the longest sub-sequence.
  -x:       Guess best frame and return the sequence with masking of the '-'
            symbols with 'X'.
  -h:       This help.

AUTHOR:
    lorenzo.cerutti\@lulix.net

EOF
exit(0);
}


my $fh;
if ($ARGV[0]) {

	open(FILE, $ARGV[0]) || die "Cannot open file $ARGV[0]\n";
	$fh = \*FILE;
}
else {

	$fh = \*STDIN;
}

my $fasta = fasta->new();
my $frames;

if ($opt_f eq '+') {
    $frames = [1,2,3];
}
elsif ($opt_f eq '-') {
    $frames = [-1,-2,-3];
}
elsif (int($opt_f) > -4 && int($opt_f) < 4) {
    $frames = [$opt_f];
}
else {
    $frames = [-3,-2,-1,1,2,3];
}

while ($fasta->to_object($fh)) {

    my @tr;

	foreach my $i (@$frames) {
		
        my $pep = $fasta->translate($i);
        $pep->set_id($fasta->get_id().".$i");
		push(@tr,$pep);
	}

    guess_best(\@tr) if ($opt_g || $opt_x);
	
    foreach (@tr) {
        
        print $_->to_fasta(); # if ($pep->get_seq());
	}
    #$fasta->delete();
}

1;

sub guess_best {
    
    my $f_list = shift;
    my $max = 0;
    my $seq = '';
    my $j;
    
    for (my $i=0; $i < scalar(@$f_list); $i++) {
        
        my $tmp_seq;
        my $s = $f_list->[$i]->get_seq();
        # count gaps
        my $g = $s =~ /\-/g;
        # count longest seq
        my $l = 0;
        foreach my $a (split /\-/,$s) {
            
            my $tmp = length($a);
            if ($tmp > $l) {
                $l = $tmp;
                $tmp_seq = $a;
            }
        }
        
        my $score = 2*$l/($l+$g);
        
        if ($score > $max) {

            $max = $score;
            $seq = $tmp_seq;
            $j = $i;
        }
    }

    my $tmp = $f_list->[$j];
    $max = int($max*100)/100;
    $tmp->set_des($tmp->get_des(). " escore=$max");
    if ($opt_x) {
        my $tmpseq = $tmp->get_seq();
        $tmpseq =~ s/-/X/g;
        $tmp->set_seq($tmpseq);
    }
    else {
        $tmp->set_seq($seq);
    }
    @$f_list = ();
    push(@$f_list, $tmp);
}
