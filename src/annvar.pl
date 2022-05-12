#!/usr/bin/env perl

use Getopt::Long;
use Data::Dumper;
use strict;

use vars qw($variants $coverage_dir $cutoff $mincvg);

GetOptions(
    "v=s"   => \$variants,
    "d=s"   => \$coverage_dir,
    "c=f"   => \$cutoff,
    "m=i"   => \$mincvg
);

$main::cutoff  = 0.66 if (!defined($cutoff));
$main::min_cvg = (defined($mincvg))? $mincvg : 10;
$main::match_score = 1;

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
         my ($variant,$score) = process_sample($db,$var,$curr,$coverage_dir);
         printf("%s,%s,%1.2f\n",$curr,$variant,$score); 
         $curr = $v[0];
         $db   = {};
    }
    
    if ($v[2] eq 'M')
    {
        my @m = split //, $v[5];
        for (my $i = 0; $i < scalar(@m); $i++)
        {
            $db->{$v[1]+$i} = [$v[2],$v[3],$m[$i]];
        }
    }
    else
    {
        $db->{$v[1]} = [$v[2],$v[3],$v[5]];
    }
}
close($fh);
exit 0;

sub read_variants
{
    my $file = shift;
    my $var  = {};
    my $name = undef;
    my $cutoff = $main::cutoff;
    
    open(my $fh, $file) || die "Cannot read VOC file.\n";
    while (<$fh>)
    {
        chomp;
     
        if (/^var\s/)
        {
            die if (!$name);

            my ($u,$num,$pos,$type,$ref,$mut,@annot) = split /\s+/;
            
            push(@{$var->{$name}{positions}{$pos}}, {
                            type => $type, 
                            num  => $num,
                            ref  => $ref,
                            mut  => $mut,
                            ann  => join(' ', @annot) });
           
            next;
        }

        if (/^name\s+(\S.*)\s*$/)
        {
            $cutoff = $main::cutoff;
            $name = $1;
            next;
        }
        if (/^cutoff\s+(\S.*)\s*$/)
        {
            $cutoff = $1;
            next;
        }
        if (/^\/\//)
        {
            $var->{$name}{cutoff} = $cutoff;
            $name = undef;
        }
    }
    close($fh);
    return $var;
}

sub process_sample
{
    my ($db,$var,$sample_name,$coverage_dir) = @_;
    my $match_score = $main::match_score;
    
    my $scores = {};
    my $cvg = [];
    if ($coverage_dir)
    {
        $cvg = read_coverage("$coverage_dir/$sample_name");
    }
    
    foreach my $voc_name (sort {$a cmp $b } keys %$var)
    {
        my $max_score = 0;
        my $score     = 0;
        my $reward    = 0;
        my $positions = $var->{$voc_name}{positions};
        my $done      = {};
        my $voc_scores= {};
        my %log;

        $log{0} = "#SAMPLE: $sample_name / $voc_name\n";
        
        foreach my $pos (sort {$a <=> $b} keys %{$positions})
        {
            my $v = $positions->{$pos};
            
            my $score_factor = ($pos >= 21563 && $pos <= 25384)? 2 : 1;
            my $reward       = $match_score * $score_factor;
            
            for (my $i = 0; $i < scalar(@{$v}); $i++)
            {
                my $v_pos = $v->[$i];
                $voc_scores->{$v_pos->{num}} = $reward;
            }
        }
        
        foreach (values %$voc_scores)
        {
            $max_score += $_;
        }

        
        foreach my $pos (sort {$a <=> $b} keys %{$positions})
        {
            my $v = $positions->{$pos};
           
            my $score_factor = ($pos >= 21563 && $pos <= 25384)? 2 : 1;
            my $reward       = $match_score * $score_factor;

            #if (!$db->{$pos})
            #{
            #    next;
            #}

            my $type  = $db->{$pos}[0];
            my $delta = $db->{$pos}[1];
            my $mut   = $db->{$pos}[2];
                
            for (my $i = 0; $i < scalar(@{$v}); $i++)
            {
                my $v_pos = $v->[$i];
                
                if ($done->{$v_pos->{num}})
                {
                    last;
                }

                my $exp = "($v_pos->{type}:$v_pos->{mut})";
                my $obs;

                if ($v_pos->{type} eq $type && $type eq 'D')
                {
                    if ($v_pos->{mut} == $delta)
                    {
                        $done->{$v_pos->{num}} = 1;
                        $score += $reward;
                        $obs = "($type:$delta)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=$reward : $score / $max_score\n";
                    }
                    elsif (!$done->{$v_pos->{num}})
                    {
                        $obs = "($type:$delta)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=0 : $score / $max_score\n";
                    }
                }
                elsif ($v_pos->{type} eq $type && $type eq 'I')
                {
                    if ($v_pos->{mut} == $delta)
                    {
                        $done->{$v_pos->{num}} = 1;
                        $score += $reward;
                        $obs = "($type:$delta)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=$reward : $score / $max_score\n";
                    }
                    elsif (!$done->{$v_pos->{num}})
                    {
                        $obs = "($type:$delta)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=0 : $score / $max_score\n";
                    }
                }
                elsif ($v_pos->{type} eq $type && $type eq 'M')
                {
                    if ($v_pos->{mut} eq $mut)
                    {
                        $done->{$v_pos->{num}} = 1;
                        $score += $reward;
                        $obs = "($type:$mut)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=$reward : $score / $max_score\n";
                    }
                    elsif (!$done->{$v_pos->{num}})
                    {
                        my $m = ($mut)? $mut : $v_pos->{ref};
                        $obs = "($type:$m)";
                        $log{$v_pos->{num}} = "#POSITIONS: $voc_name,$pos,$cvg->[$pos]".",obs=$obs,exp=$exp,$v_pos->{ann},SCORE=0 : $score / $max_score\n";
                    }
                }
            }
            
            for (my $i = 0; $i < scalar(@{$v}); $i++)
            {
                my $v_pos = $v->[$i];
                
                if ($done->{$v_pos->{num}})
                {
                    next;
                }
                else 
                {
                    my $message = ($cvg->[$pos] < $main::min_cvg)? 
                                    "NO_COVERAGE" : "MISSING_MUTATION";
                    my $exp = "($v_pos->{type}:$v_pos->{mut})";
                    my $obs = ($cvg->[$pos] < $main::min_cvg)? "n.a." : 
                              (defined($mut))? "($type:$mut)" : "n.a.";
            
                    $log{$v_pos->{num}} = 
                        "#POSITIONS: $voc_name,$pos,$cvg->[$pos]" .
                        ",obs=$obs,exp=$exp,$v_pos->{ann},$message\n";
                } 
            }
        }
        
        $scores->{$voc_name} = $score/$max_score;
        $score = sprintf("%1.2f", $score/$max_score);
        
        foreach (sort {$a <=> $b} keys %log)
        {
            print STDERR $log{$_};
        }
        print STDERR "#SCORE: $voc_name,$score\n";
        print STDERR "#//\n";
    }
    
    my $attribution = 'not_attributed';
    my $fscore      = 0;
    foreach my $name (keys %$scores)
    {
        if ($scores->{$name} > $fscore && $scores->{$name} >= $var->{$name}{cutoff})
        {
            $attribution = $name;
            $fscore      = $scores->{$name};
        }
    }
    return ($attribution,$fscore);
}

sub read_coverage
{
    my $file = shift;
    my $filename = `ls $file.cvg`; chomp($filename);
    my $cvg = [];
    open(my $cf, $filename) || die "Cannot open '$filename'\n";
    while (<$cf>)
    {
        my ($pos,$c) = /^\S+\s+(\d+)\s+(\d+)/;
        $cvg->[$pos] = $c;
    }
    return $cvg;
}


__END__


