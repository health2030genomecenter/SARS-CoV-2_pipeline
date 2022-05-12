#!/usr/bin/env perl

use constant REFERENCE => "NC_045512";
use constant INTERVAL  => 16;

our %f = (
    'fwd' => '0x42',
    'rev' => '0x82'
);

use Data::Dumper;
use strict;

my $bam = $ARGV[0];
my $pos = $ARGV[1];
my $ref = REFERENCE;

my $h = {};
my $c = {};
my $x = parse_bam('fwd',$h,$pos,$c);
$x   .= parse_bam('rev',$h,$pos,$c);

printf("%7s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\n",
    "Var",">-d","<-d","-d",">+d","<+d","+d");

foreach my $var (sort {$a cmp $b} keys %$h)
{
    printf("%5i:%s\t%5i\t%5i\t%5i\t%5i\t%5i\t%5i\n", 
        $pos,$var,
        $h->{$var}{fwd}{wo_dup},
        $h->{$var}{rev}{wo_dup},
        $h->{$var}{tot}{wo_dup},
        $h->{$var}{fwd}{w_dup},
        $h->{$var}{rev}{w_dup},
        $h->{$var}{tot}{w_dup});
}

my $max = 1;
foreach my $p (sort {$a <=> $b} keys %$c)
{
    $max = $c->{$p} if ($c->{$p} > $max);
}

print "\nCoverage:\n";
for (my $p = $pos - INTERVAL; $p <= $pos + INTERVAL; $p++)
{
    my $cvg = $c->{$p} || 0;
    printf("%8i : %8i",$p,$cvg);
    ($p == $pos)? print " $x|" : print "   |";
    my $x = int($cvg/$max * 25);
    print "=" x $x;
    print "\n";
}

exit 0;

sub parse_bam
{
    my ($f,$h,$pos,$cvg) = @_;
    my $filter = $f{$f};
    my $cmd = "samtools view $bam $ref:$pos-$pos -f $filter";
    my %support;
    my %highest;

    open(my $fh, "$cmd |") || die "Cannot execute $cmd\n";
    while (<$fh>)
    {
        chomp; 
        my @v      = split/\s+/; 
        my $s      = $v[3]; 
        my $n      = 0;
        my @x      = split //,$v[9];
        my $k      = $s.":".$v[5];

        LOOP: while ($v[5] =~ s/^(\d+)([MDI])//)
        {
            if ($2 eq 'M')
            {
                for (my $i = 0; $i < $1; $i++)
                {
                    if ($s == $pos)
                    {
                        $support{$x[$n]}{$k}++;
                        #last LOOP;
                    }

                    if ($s >= $pos-INTERVAL && $s <= $pos+INTERVAL)
                    {
                        $cvg->{$s}++;
                    }
                    
                    $n++;
                    $s++;
                }
            }
            elsif ($2 eq 'D')
            {
                for (my $i = 0; $i < $1; $i++)
                {
                    if ($s == $pos)
                    {
                        $support{'-'}{$k}++;
                        #last LOOP;
                    }
                    $s++;
                }
            }
            elsif ($2 eq 'I')
            {
                for (my $i = 0; $i < $1; $i++)
                {
                    if ($s == $pos)
                    {
                        $support{'.'}{$k}++;
                        #last LOOP;
                    }
                    $n++;
                }
            }
            last LOOP if ($s == $pos+INTERVAL);
        }
    }

    foreach my $var (sort {$a cmp $b} keys %support)
    {
        my $without_dup = scalar(keys %{$support{$var}});
        my $with_dup    = 0;
        foreach (keys %{$support{$var}})
        {
            $with_dup += $support{$var}{$_};
        }
        $h->{$var}{$f}{'wo_dup'}  += $without_dup;
        $h->{$var}{$f}{'w_dup'}   += $with_dup;
        $h->{$var}{tot}{'wo_dup'} += $without_dup;
        $h->{$var}{tot}{'w_dup'}  += $with_dup;
    }
    
    my $max = 0;
    my $b   = '*';
    foreach my $var (sort {$a cmp $b} keys %support)
    {
        if ($h->{$var}{tot}{'wo_dup'} > $max)
        {
            $max = $h->{$var}{tot}{'wo_dup'};
            $b   = $var;
        }
    }
    return $b;
}
