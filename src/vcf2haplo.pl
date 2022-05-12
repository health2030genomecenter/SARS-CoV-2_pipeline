#!/usr/bin/env perl

use constant FILTER => 650;

use Getopt::Long;
use Data::Dumper;
use strict;

my %opt;
GetOptions( \%opt, "filter=i" );

my $filter = $opt{filter} || FILTER;
warn("Filter=$filter\n");

my %db; # keep informatin needed for report
my @ht; # observed haplotypes (each haplotype is represented by string)

my $fh = \*STDIN;

while (<$fh>)
{
    next if (!/^NC_045512/);

    my @v = split /\s+/;
    my ($pos,$ref,$var,$p) = ($v[1],$v[3],$v[4],$v[5]);
    
    next if ($filter > $p); 
    
    my $vars = []; # store reference and observed variants
    foreach my $v (split /,/, $var)
    {
        push(@$vars,$v);
    }
    unshift(@$vars,$ref); 
    
    $db{$pos}{var} = $vars;

    my @w = split /\;/, $v[7];
    my @z = split /\:/, $v[9];
    $w[0] =~ s/AB=//;
    
    my $freqs = []; # store reference and observed variants freqs
    my $sum   = 0;
    foreach my $f (split /,/, $w[0])
    {
        $f = 1 if ($f == 0);
        push(@$freqs,$f);
        $sum += $f;
    }
    unshift(@$freqs,1-$sum);

    # Generate the haplotype string based on observed variants
    my $i = 0;
    foreach (split /[\|\/]/, $z[0])
    {
        if ($_ eq '.')
        {
            $_ = 1;
        }
        @ht[$i] .= $_;
        $i++;
    }
    
    $db{$pos}{frq} = $freqs;
}
close($fh);

# Collapse haplotypes
my %ht = map { $_ => 1 } @ht;

my @uniq; # Counts of haplotypes having the variant at a given pos
foreach my $ht (keys %ht)
{
    my $i = 0;
    foreach (split //,$ht)
    {
        $uniq[$i]{$_}++;
        $i++;
    }
}
my @nh; # Get observed haplotypes for a given position
for (my $i = 0; $i < scalar(@uniq); $i++)
{
    $nh[$i] = scalar(keys %{$uniq[$i]});
}

# Report
my %hf;
my %ho;
my $h = 0;
foreach my $ht (keys %ht)
{
    my @ht  = split //, $ht;
    my @pos = sort {$a <=> $b} keys %db; 

    if (scalar(@ht) != scalar(@pos))
    {
        die "Inconsistent genotype/variant positions.\n";
    }
    
    my $min = 0;
    my $n   = 0;
    my $hu  = '';
    for (my $i = 0; $i < scalar(@ht); $i++)
    {
        my $gt  = $ht[$i];
        my $mut = $db{$pos[$i]}{var};
        my $frq = $db{$pos[$i]}{frq};
        my $dis = ($nh[$i] > 1)? '*' : ' ';
        $hu .= $gt if ($dis eq '*');

        $ho{$ht} .= sprintf("%3s %1s %6i\t%-4.2f%\t%s -> %s\n",$gt,$dis,$pos[$i],100*$frq->[$gt],$mut->[0],$mut->[$gt]);
    
        # black magic estimation of haplotype freqs (but probably wrong)
        if ($frq->[$gt] < 1 && $mut->[0] && $nh[$i] > 1)
        {
            $min += $frq->[$gt];
            $n += $uniq[$i]{$gt};
        }
    }
    
    my $m = ($n)? $min/$n : 1;
    $h++;
    print ">Haplotype * $h: ",sprintf("%1.2f",$m)," [$hu]\n";
    print $ho{$ht};
    print "//\n";
}
