#!/usr/bin/env perl

use constant REFERENCE => "NC_045512";

use strict;

my $bam = $ARGV[0];
my $pos = $ARGV[1];
my $mut = $ARGV[2];
my $ref = REFERENCE;

my ($f_support,$fd_support,$f_best,$f_bestscore,$f_t) = parse_bam('0x42');
my ($r_support,$rd_support,$r_best,$r_bestscore,$r_t) = parse_bam('0x82');

printf("Fow(withDup) = %1.2f \t Fow(withoutDup) = %1.2f \t Best = %s \t Best_score = %1.2f \n",$f_support,$fd_support, $f_best, $f_bestscore,$f_t);
printf("Rev(withDup) = %1.2f \t Rev(withoutDup) = %1.2f \t Best = %s \t Best_score = %1.2f \n",$r_support,$rd_support, $r_best, $r_bestscore,$r_t);

my $sig       = 1/(1+exp(-100*($f_t+$r_t-100)));
my $g_support = $f_support + $fd_support; 
$g_support   += $r_support + $rd_support; 
$g_support   /= 4;
printf("Support = %1.2f% (%i)\n",100*$sig*$g_support,$f_t+$r_t);

exit 0;

sub parse_bam
{
    my $filter = shift;
    my $cmd = "samtools view $bam $ref:$pos-$pos -f $filter";
    my %support;
    my %highest;
    my ($o,$t) = (0,0);
    open(my $fh, "$cmd |") || die "Cannot execute $cmd\n";
    while (<$fh>)
    {
        chomp; 
        my @v      = split/\s+/; 
        my ($n,$m) = $v[5] =~ /^(\d+)M(\d+\D)/; 
        my $s      = $v[3]+$n; 
        my $str    = $s."\t".$m."\t".$v[3]."\t".$v[5];

        if ($s == $pos && $m eq $mut)
        {
            $support{o}{$str}++;
            $o++;
        }
        else
        {
            $support{n}{$str}++;
        }

        $highest{"$s $m"}++;
        $support{t}{$str}++;
        $t++;
    }
   
    foreach (keys %{$support{o}})
    {
        print STDERR "SUPPORTING_OBS($filter): $support{o}{$_} \t $_\n";
    }
    foreach (keys %{$support{n}})
    {
        print STDERR "NOT_SUPPORTING_OBS($filter): $support{n}{$_} \t $_\n";
    }
    
    my ($best, $best_score);
    foreach (sort {$highest{$b} <=> $highest{$a}} keys %highest)
    {
        $best = $_;
        $best_score = $highest{$_}/$t;
        last;
    }

    my $full_s = ($o)? $o/$t : 0;
    my $dup_s  = ($o)? scalar(keys %{$support{o}})/scalar(keys %{$support{t}}) : 0;

    return ($full_s, $dup_s, $best, $best_score,$t);
}
