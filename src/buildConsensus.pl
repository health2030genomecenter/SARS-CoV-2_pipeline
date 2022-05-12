#!/usr/bin/env perl

use constant MIN_CVG     => 10;
use constant MIN_MUT_CVG => 10;
use constant WIDTH       => 5;

use Data::Dumper;
use Getopt::Long;
use strict;

my $ref;
my $cov;
my $min;
my $min_mut;
my ($s,$e);

GetOptions(
    "f=s" => \$ref,
    "c=s" => \$cov,
    "m=i" => \$min,
    "n=i" => \$min_mut,
    "s=i" => \$s,
    "e=i" => \$e
);

my ($id,$refseq) = read_ref($ref);
my @ref  = split //, $refseq;
my $cov  = read_cov($cov);
$min     = MIN_CVG if (!defined($min));
$min_mut = MIN_MUT_CVG if (!defined($min_mut));

# process vcf
my %vcf = {};
while (<>)
{
    chomp;
    if (/^$id\s/)
    {
        my @v = split /\s+/;
        my $pos = $v[1] - 1;
        $vcf{$pos}{ref} = $v[3];
        $vcf{$pos}{var} = $v[4];
        $vcf{$pos}{mod} = (length($v[3]) > length($v[4]))? 'D' :
                           (length($v[3]) < length($v[4]))? 'I' : 'M';
    }
}

# build consensus
my $cons;
for (my $i = 0; $i < scalar(@$cov); $i++)
{
    if ($s && $e)
    {
        next if ($i < $s-1);
        next if ($i > $e-1);
    }

    if ($vcf{$i} && $cov->[$i] >= $min_mut)
    {
        if ($vcf{$i}{mod} eq 'M')
        {
            my @var = split //, $vcf{$i}{var};
            my @rm  = split //, $vcf{$i}{ref};
            my $n=0;
            foreach (@var)
            {
                $cons .= $_;
                $n++;
            }
            $i += $n-1;
        }
        elsif ($vcf{$i}{mod} eq 'I')
        {
            my @ins = split //, $vcf{$i}{var};
            my @rm  = split //, $vcf{$i}{ref};
            foreach (@ins)
            {
                $cons .= $_;
            }

            $i += scalar(@rm)-1;
        }
        elsif ($vcf{$i}{mod} eq 'D')
        {
            my @del = split //, $vcf{$i}{var};
            my @rm  = split //, $vcf{$i}{ref};
            foreach (@del)
            {
                $cons .= $_;
            }

            $i += scalar(@rm)-1;
        }
    }
    elsif ($vcf{$i})
    {
        $cons .= 'N';
    }
    else
    {
        $cons .= ($cov->[$i] < $min)? 'N' : $ref[$i];
    }
    #print "$i -> ",substr($cons,-20),"\t",$vcf{$i}{ref},"->",$vcf{$i}{var},":",$vcf{$i}{mod},"\n";
}

# print consensus
print ">$id\n";
print "$1\n" while ($cons =~ /(.{1,60})/g);

exit 0;

sub read_cov
{
    my $cov = [];
    open(my $fh, $_[0]) || die "Cannot open coverage -c $_[0]\n";
    while (<$fh>)
    {
        chomp;
        my ($unused,$i,$c) = split /\s+/;
        $cov->[$i] = $c;
    }
    close($fh);
    return $cov;
}

sub read_ref
{
    my $id;
    my $seq;

    open(my $fh, $_[0]) || die "Cannot open reference -f $_[0]\n";
    while (<$fh>)
    {
        if (/^>(\S+)/)
        {
            $id = $1;
            next;
        }

        chomp;
        $seq .= $_;
    }
    close($fh);
    return($id,$seq);
}


