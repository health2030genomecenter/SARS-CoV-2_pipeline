#!/usr/bin/env perl

use Getopt::Long;
use Data::Dumper;
use strict;

use vars qw($variants $coverage_dir $cutoff);

GetOptions(
    "v=s"   => \$variants,
    "d=s"   => \$coverage_dir,
    "c=f"   => \$cutoff
);

$cutoff = 0.66 if (!defined($cutoff));

# read variants description
my $var = read_variants($variants);

my $curr;
my $db = {};
my $fh = \*STDIN;

while (<$fh>)
{
    chomp;
    my @v = split /,/;

    if (!$curr)
    {
        $curr = $v[0];
    }

    if ($curr ne $v[0] || eof($fh))
    {
        my ($variant,$score) = process_sample($db,$var,$cutoff,$curr,$coverage_dir);
        printf("%s,%s,%1.2f\n",$curr,$variant,$score); 
        $curr = $v[0];
        $db   = {};
    }
    
    $db->{$v[1]} = [$v[2],$v[3],$v[5]];
}
close($fh);
exit 0;

sub process_sample
{
    my ($db,$var,$cutoff,$sample_name,$coverage_dir) = @_;
    my $scores = {};
    
    my $cvg = [];
    if ($coverage_dir)
    {
        $cvg = read_coverage("$coverage_dir/$sample_name");
    }

    foreach my $name (sort {$a cmp $b } keys %$var)
    {
        print STDERR "#SAMPLE: $sample_name\n";
        my $gscore     = 0;
        my $v = $var->{$name}{pos};
        my $t = $var->{$name}{tot_mut};

        foreach my $pos (sort {$a <=> $b} keys %$v)
        {
            my $a = $var->{$name}{annot}{$pos};
            if ($db->{$pos})
            {
                my $type  = $db->{$pos}[0];
                my $delta = $db->{$pos}[1];
                my $mut   = $db->{$pos}[2];
            
                my $score     = 1;
                my $tot_score = 3;
                my $str;
                if ($v->{$pos}{type} eq $type && $type eq 'D')
                {
                    if ($v->{$pos}{mut} == $delta)
                    {
                        $score += 1;
                    }
                    $score += 1;
                    $str    = "($type,$delta)";
                }
                elsif ($v->{$pos}{type} eq $type && $type eq 'I')
                {
                    if ($v->{$pos}{mut} == $delta)
                    {
                        $score += 1;
                    }
                    $score += 1;
                    $str    = "($type,$delta)";
                }
                elsif ($v->{$pos}{type} eq $type && $type eq 'M')
                {

                    if (index($mut,$v->{$pos}{mut}) > -1)
                    {
                        $score += 1;
                    }
                    $score += 1;
                    $str    = "($type,$mut)";
                }
                
                $score      /= $tot_score;
                $gscore     += $score;
                
                
                print STDERR "#POSITIONS: $name,$pos,$cvg->[$pos-1]",
                        ",obs=$str",
                        ",exp=(",$v->{$pos}{type},",",$v->{$pos}{mut},
                        "),$a\n";
            }
            else
            {
                my $str   = "na";
                print STDERR "#POSITIONS: $name,$pos,$cvg->[$pos-1]",
                        ",obs=$str",
                        ",exp=(",$v->{$pos}{type},",",$v->{$pos}{mut},
                        "),MISSING\n";
            }
        }
        if ($gscore > 0)
        {
            printf(STDERR "#SCORE: %s,%1.2f\n",$name,$gscore/$t);
            $scores->{$name} = $gscore/$t;
        }
        else
        {
            printf(STDERR "#SCORE: %s,%1.2f\n",$name,0);
            $scores->{$name} = 0;
        }
        print STDERR "#//\n";
    }
    
    my $attribution = 'unknown';
    my $fscore      = 0;
    foreach my $name (keys %$scores)
    {
        if ($scores->{$name} > $fscore && $scores->{$name} >= $cutoff)
        {
            $attribution = $name;
            $fscore      = $scores->{$name};
        }
    }
    return ($attribution,$fscore);
}

sub read_variants
{
    my $file = shift;
    my $var  = {};
    my $name = undef;
    my $sum  = {};

    open(my $fh, $file) || die "Cannot read variant description file.\n";
    while (<$fh>)
    {
        chomp;
        if (/^var\s/)
        {
            die if (!$name);
            my ($u,$num,$pos,$type,$mut,@annot) = split /\s+/;
            $sum->{$name}{$num}++;
            $var->{$name}{tot_mut}++ if ($sum->{$name}{$num} == 1);
            $var->{$name}{pos}{$pos} = {type => $type, mut => $mut};
            $var->{$name}{annot}{$pos} = join ' ', @annot;
            next;
        }
        if (/^name\s+(\S.*)\s*$/)
        {
            $name = $1;
            next;
        }
        if (/^\/\//)
        {
            $name = undef;
        }
    }
    close($fh);
    return $var;
}

sub read_coverage
{
    my $file = shift;
    my $filename = `ls $file*`; chomp($filename);
    my $cvg = [];
    open(my $cf, $filename) || die "Cannot open $filename\n";
    while (<$cf>)
    {
        my ($pos,$c) = /^(\d+)\s+(\d+)/;
        $cvg->[$pos] = $c;
    }
    return $cvg;
}
