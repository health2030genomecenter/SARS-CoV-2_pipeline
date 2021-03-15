package BioSeqFormat;
$RELEASE=0.7.4;


#-------------------------------------------------------------#	
###############################################################	
###############################################################
#
# This module is developed and maintained by Lorenzo CERUTTI
# (lorenzo.cerutti@lulix.net).
#
# You can use it and modify it, just leave this comment and 
# author information. ... be nice with the author :)
#
###############################################################
###############################################################
#-------------------------------------------------------------#	


#-------------------------------------------------------------#
# RELEASE 0.7.4
#	- change checksum calculation based on md5_base64
#
# RELEASE 0.7.3
#	- modify to_object methof of fasta to read from a pipeline
#	  and add function peek() to read a char from a buffer without
#	  removing it from the buffer.
#
# RELEASE 0.7.2
#	- add read() which is redundant to to_object().
#	- add pos(), delete_pos() and identity() methods.
#
# RELEASE 0.7.1
#	- add window_shuffle method to randomize sequences.
# 	- add pos() and delete_pos() methods
#
# RELEASE 0.7
#	- add PREDATOR secondary structure prediction method: predator()
#	- modify JNET prediction (jnet()) to behave as predator(): returns
#	  an object.
#	- add ungap() and gap() methods to delete and reset gaps in a 
#	  sequence.
#
# RELEASE 0.6
#	- complement(): now modify directly the object and returns void
#	- reverse(): now modify directly the object and returns void
#	- translate(): now returns an object with the tranlated sequence.
#
# RELEASE 0.5 (23/09/2002)
#	- to_object() in fasta class modified to accept fasta IDs with
#	  a \s* at the end.
#
# RELEASE 0.4
#	- Modify the to_object fasta (some problems with some files).
#-------------------------------------------------------------#	

use Carp;
use GenObj;
use SeqCodes;
use IO::Handle;
#use POSIX ":stdio_h";
#use Math::Random;
use Digest::MD5  qw(md5_hex);
use strict;

@BioSeqFormat::ISA = qw(GenObj);
@BioSeqFormat::_ATTRIBUTES_ = qw(

	id:val
	ac:list 
	des:val
	seq:val 
	qual:val
    type:val 
	length:val 
	checksum:val 
	ft:list
	dr:list
    oc:val
	ox:val
    os:val
    zv:val
);

###############################################################
# Internal functions
###############################################################

sub array_fisher_yates_shuffle 
{
	my $array = shift; 
	my $i; 
	for ($i = @$array; --$i; ) 
	{ 
		my $j = int rand ($i+1); 
		next if $i == $j; 
		@$array[$i,$j] = @$array[$j,$i]; 
	}
}

#-------------------------------------------------------------#
# This function reproduces the peek function of C++
#-------------------------------------------------------------#	
sub peek($) {

	my $c  = getc($_[0]);
	IO::Handle::ungetc($_[0],ord($c));
	return $c;
}


###############################################################
# General methods
###############################################################


#-------------------------------------------------------------#
#
# guess_type()
#
# DESCRIPTION:
#	To guess the type of the sequence (P|N|?)
#
# USAGE:
#	BioSeqFormat::guess_type($string);
#
# ARGUMENTS:
#	$string: Contains the sequence
#
# RETURN:
# 	char: 'P' for protein, 'N' for nucleotides
#
#-------------------------------------------------------------#
sub guess_type($) {
	
	my $seq = shift;
	
	if ($seq =~ /[QWERYIPSDFHKLVNMX]/g) {
		
		return 'P';
	}
	elsif ($seq =~ /[ACGTN]/g) {
		
		return 'N';
	}
	else {
		return '?';
	}
}

#-------------------------------------------------------------#
#
# checksum()
#
# DESCRIPTION:
#	Return a simplified checksum
#
# USAGE:
#	BioSeqFormat::checksum($string);
#   $object->checksum();
#
# ARGUMENTS:
#	$string: Contains the sequence
#
# RETURN:
# 	The checksum value
#
#-------------------------------------------------------------#
sub checksum($) {
	
	if (ref($_[0])) {
	
		return md5_hex($_[0]->get_seq());
		#return (unpack("%64C*",$_[0]->get_id() . $_[0]->get_seq()))
	}
	else {
		
		return md5_hex($_[0]);
		#return (unpack("%64C*",$_[0]));
	}
}
	
#-------------------------------------------------------------#
#
# aa3to1
#
# Description:
#	To convert a 3 letter code aa in a 1 letter code
#
# Usage:
#	BioSeqFormat::aa3to1($string);
#
# Arguments:
#	$string: aa in 3 letter code
#
# Return:
#	string of the 1 letter code aa
#
#-------------------------------------------------------------#
sub aa3to1() {

	return $SeqCodes::to_1_letter_code_aa{uc($_[0])};
}

#-------------------------------------------------------------#
#
# aa1to3
#
# Description:
#	To convert a 1 letter code aa in a 3 letter code
#
# Usage:
#	BioSeqFormat::aa1to3($string);
#
# Arguments:
#	$string: aa in 1 letter code
#
# Return:
#	string of the 3 letter code aa
#
#-------------------------------------------------------------#
sub aa1to3() {

	return $SeqCodes::to_3_letter_code_aa{uc($_[0])};
}


#-------------------------------------------------------------#
#
# set_seed()
#	Set a seed for the random number generator	
#
# USAGE:
#	$fasta->set_seed($seed);
#
# ARGUMENTS:
#	$seed:   seed for the random generator
#
# RETURNS:
#	void
#
#-------------------------------------------------------------#
sub BioSeqFormat::set_seed($) {

	srand($_[1]);
	#Math::Random::random_set_seed_from_phrase($_[1]);
}

#-------------------------------------------------------------#
#
# generate_rnd_seq()
#	Generate a random sequence in the passed object.	
#
# USAGE:
#	$fasta->generate_rnd_seq($freq_hash,$length,$seed);
#
# ARGUMENTS:
#	$freq_hash: hash containing the frequencies
#	$length:    integer representing the length of the sequence
#	            to be generated
#	$seed:      seed for the random generator
#
# RETURNS:
#	void
#
#-------------------------------------------------------------#
sub BioSeqFormat::generate_rnd_seq(@) {

	my ($self,$freq_hash,$len,$seed) = @_;
	my $rnd_seq;
	my ($max,@sum_freq,@alphabet);
	
	foreach (keys %{$freq_hash}) {
		
		$max += $freq_hash->{$_};
		push(@sum_freq,$max);
		push(@alphabet,$_);
	}

	$self->set_seed($seed) if (defined($seed));
	
	my ($i,$j);
	#FIXME
	my @r = ();#Math::Random::random_uniform($len,0,$max);
	
	for ($i = 0; $i < $len; $i++) {
		
		LOOP:for ($j = 0; $j < scalar(@alphabet); $j++) {
			
			if ($r[$i] < $sum_freq[$j]) {
				
				$rnd_seq .= $alphabet[$j];
				last LOOP;
			}
		}
	}

	$self->set_seq($rnd_seq);
}


#-------------------------------------------------------------#
#
# translate() 
#
# Description:
#	To translate a DNA sequence in protein
# Usage:
#	$new_object = $object->translate(int,bool);
#
# Arguments:
#	int: frame [1,2,3] plus strand; [-1,-2,-3] minus strand.
#   bool: if true use 'fs' to indicate frameshifts.
#
# Return:
#	A translated fasta object.
#
#-------------------------------------------------------------#
sub translate() {

   my ($self,$frame,$informFS) = @_;
   my $i;

   return undef if ($frame < -3 || $frame > 3 || !defined($frame));

   my $clone = $self->clone();

   if ($frame < 0) {

      $clone->reverse();
      $clone->complement();
      $frame *= -1;
   }	

   my @seq = split //,$clone->get_seq();
   my $translation;
   $frame--;

   for ($i=$frame;$i<(scalar(@seq)-2);$i+=3) 
   { 
      my $lower = (ord($seq[$i]) > 96 && ord($seq[$i+1]) > 96 && ord($seq[$i+2]) > 96);
      my $triplet = uc($seq[$i].$seq[$i+1].$seq[$i+2]);
      my $aa;
      if ($informFS)
      {
         if (0 == length($triplet) % 3 && $triplet !~ /[^A-Za-z]/)
         {
            $aa = ($SeqCodes::codon{$triplet})? $SeqCodes::codon{$triplet} : 'X';
         }
         else
         { 
            $aa = '!';
         }
      }
      else
      {
         if ($triplet =~ /\-/)
         { 
            $aa = '-';
         }
         else
         {
            $aa = ($SeqCodes::codon{$triplet})? $SeqCodes::codon{$triplet} : 'X';
         }
      }
      $translation .= ($lower)? lc($aa) : $aa;
   }
   $clone->set_seq($translation);
   return $clone;
}


sub translate_codon() {

   my ($self,$informFS) = @_;
   my $i;

   my $clone = $self->clone();

   my $seq   = $self->get_seq();
   my $seqUC = $seq; $seqUC =~ s/[a-z\-]//g;   
   my $aa;

   if (length($seq) == 3)
   {
      $aa = ($SeqCodes::codon{$seqUC})? $SeqCodes::codon{$seqUC} : 'X';
   }
   else
   {
      my $tmp = uc($seq);
      while ($tmp =~ s/^(.{3})//)
      {
         $aa .= ($SeqCodes::codon{$1})? $SeqCodes::codon{$1} : 'X';
      }
   }
   
   if ($seq =~ /\-/)
   {
      $aa = ($seq eq '---')? '-' : 
            ($informFS    )? '!' : '-';
   }
   elsif (length($seq) > 3)
   {
      $aa = ($informFS)? '!' : $aa;
   }
   $clone->set_seq($aa);
   return $clone;
}
#-------------------------------------------------------------#
#
# aa2dna() 
#
# Description:
#	To translate a DNA sequence in protein
#
# Usage:
#	$new_object = $object->translate(PACKED=>(0|1));
#
# Arguments:
#	PACKED => (0|1)
#	                where 0 => max 1 IUPAC per codon (random choose between multiple codons)
#	                      1 => 1 codon per aa (codon can contain multiple IUPAC symbols)
#
# Return:
#	A dna sequence object.
#
#-------------------------------------------------------------#
sub aa2dna() 
{
    my $self = shift;
    my %opt  = @_;
    my $translator = ($opt{PACKED})? \%SeqCodes::aa2dna_packed : \%SeqCodes::aa2dna;
    my $same_codon = ($opt{SAME_CODON}) || 0;

	my $clone = $self->clone();
	my @seq = split //,$clone->get_seq();
	my $translation;
    
    if ($opt{PACKED})
    {
        foreach (@seq)
        {
            $translation .= $translator->{$_};
        }
    }
    else
    {
        foreach (@seq)
        {
            my $n = scalar(@{$translator->{$_}});
            my $i = ($same_codon)? 0 : int($n*rand(1));
            $translation .= $translator->{$_}->[$i];
        }
    }
	$clone->set_seq($translation);
	return $clone;
}

#-------------------------------------------------------------#
#
# jnet
#
# Description:
#	To predict secondary structure (REQUIRES JNET!!)
#
# Usage:
#	$new_object = $object->jnet();
#
# Return:
#	Object containing secondary structure prediction
#
#-------------------------------------------------------------#
sub jnet() {

	my $self = shift;
	my $jnet = `which jnet` || die "ERROR: jnet prediction requires jnet\n";;
	my $tmp  = GenObj::TMP . "jnet$$.tmp";
	
	chomp($jnet);
	
	my $ungap_seq = $self->clone();
	my $v = $ungap_seq->ungap();
	open(TMP,">$tmp") || die "ERROR: Cannot create tmp file\n";	
	print TMP $ungap_seq->to_fasta();
	close(TMP);
	
	my $sstr;
	open(PRED,"$jnet -c $tmp 2>/dev/null |") || die "ERROR: Cannot execute $jnet\n";

	while(<PRED>) {

		if (/jnetpred:(\S+)\,\n/) {
			
        	$sstr = $1;
			$sstr =~ s/\,//g;
			$sstr =~ s/\-/$SeqCodes::structure_code{'coil'}/g;
			$sstr =~ s/H/$SeqCodes::structure_code{'helix'}/g;
			$sstr =~ s/E/$SeqCodes::structure_code{'extend'}/g;
		}
	}	
	
	close(PRED);
	unlink($tmp);
	my $pred = $self->clone();
	$pred->set_seq($sstr);
	$pred->gap($v);
	return (defined($sstr))?$pred:undef;
}

#-------------------------------------------------------------#
#
# predator()
#
# Description:
#	To predict secondary structure (REQUIRES PREDATOR!!)
#
# Usage:
#	$new_object = $object->predator();
#
# Return:
#	Object containing secondary structure.
#
#-------------------------------------------------------------#
sub predator() {

	my $self     = shift;
	my $predator = `which predator` || die "ERROR: predator prediction requires predator\n";
	my $tmp      = GenObj::TMP . "predator$$.tmp";
	
	chomp($predator);
	
	my $ungap_seq = $self->clone();
	my $v = $ungap_seq->ungap();
	open(TMP,">$tmp") || die "ERROR: Cannot create tmp file\n";	
	print TMP $ungap_seq->to_fasta();
	close(TMP);
	
	my $sstr;
	open(PRED,"$predator $tmp |") || die "ERROR: Cannot execute predator\n";
	
	while(<PRED>) {
		
		chomp;

		if (/^\s+([_HE]+)/) {
			
			my $tmp = $1;
			$tmp =~ s/_/$SeqCodes::structure_code{'coil'}/g;
			$tmp =~ s/H/$SeqCodes::structure_code{'helix'}/g;
			$tmp =~ s/E/$SeqCodes::structure_code{'extend'}/g;
			$sstr .= $tmp;
		}
	}	
	close(PRED);
	unlink($tmp);
	my $pred = $self->clone();
	$pred->set_seq($sstr);
	$pred->gap($v);
	return (defined($sstr))?$pred:undef;
}

#-------------------------------------------------------------#
#
# ungap()
#
# Description:
#	Clear gaps of a sequence ('-' chars) and return a ref to a 
#	vector of the coordinates of the deleted gaps.
#
# Usage:
#	$vect = $object->ungap();
#
# Return:
#	Reference to vector containing the coordinates of the deleted gaps
#
#-------------------------------------------------------------#
sub ungap() {

	my $self = shift;
	my @v;
	my @seq = split('',$self->get_seq());
	my $len = scalar(@seq);
	my $i;
	my $ungap_seq;

	for ($i = 0; $i < $len; $i++) {
		
		if ($seq[$i] eq '-' || $seq[$i] eq '.') {

			push(@v,$i);
		}
		else {

			$ungap_seq .= $seq[$i];
		}
	}
	
	$self->set_seq($ungap_seq);
	return(\@v);
}

#-------------------------------------------------------------#
#
# gap()
#
# Description:
#	Add gaps to a sequence ('-' chars). Normally this function is used
#	together with the 'ungap' function.
#
# Usage:
#	$object->gap(\@vect);
#
# Arguments:
#	\@vec: reference list containing the coordinates of the gaps
#
# Return:
#	void
#
#-------------------------------------------------------------#
sub gap() {

	my ($self,$v) = @_;
	my $len = scalar(@$v);
	my @seq = split //, $self->get_seq();
	my $i;

	for ($i = 0; $i < $len; $i++) {
		
		my @tmp = @seq;
		my $j;
		$seq[$v->[$i]] = '-';
		
		for ($j = $v->[$i]; $j < scalar(@tmp); $j++) {

			$seq[$j+1] = $tmp[$j];
		}
	}
	
	$self->set_seq(join('',@seq));
}

#-------------------------------------------------------------#
# 
# pos()
#
# Description:
#	Returns symbol at the position specified of the sequence.
#
# Usage:
#	my $char = $object->pos($int,$char);
#
# Arguments:
#	$int:  the position in the sequence (1..len)
#	$char: the symbol to set at the specified position
#
# Returns:
#	$char: the symbol at the specified position
#
#-------------------------------------------------------------#
sub pos() {

	my ($self,$pos,$symbol) = @_;

	if ($pos < 1 || $pos > $self->get_length()) {

		return undef;
	}

	if ($symbol) {
		
		my $seq = $self->get_seq();
		substr($seq,$pos-1,1) = $symbol;
		$self->set_seq($seq);
	}
	
	return (substr($self->get_seq(),$pos-1,1));
}

#-------------------------------------------------------------#
# 
# delete_pos()
#
# Description:
#	Delete specified positions in the sequence.
#
# Usage:
#	$object->delete_pos(@list);
#
# Arguments:
#	@list: the positions in the sequence (1..len) to delete
#
# Returns:
#	void
#
#-------------------------------------------------------------#
sub delete_pos() {

	my ($self,@pos) = @_;
	my $new_seq;
	my $len = $self->get_length();

	LOOP: for (my $i = 1; $i <= $len; $i++) {
		
		foreach (@pos) {

			next LOOP if ($_ == $i);
		}
		
		$new_seq .= $self->pos($i);
	}

	$self->set_seq($new_seq);
}

#-------------------------------------------------------------#
# 
# identity()
#
# Description:
#	Return identity % between two aligned sequences.
#
# Usage:
#	float $p_ident = $object1->identity($object2);
#
# Arguments:
#	$object2: a sequence object
#
# Returns:
#	float: % identity between the 2 object sequences (0-1) or 
#	       UNDEF if sequences are not aligned (different length).
#
#-------------------------------------------------------------#
sub identity($) : method {

	my ($seq1,$seq2) = @_;
	my $length = $seq1->get_length();

	if ($length != $seq2->get_length()) {

		print STDERR "Could not mesure % identity between sequences. Different length!\n";
		return undef;
	}

	my $identity = 0;
	my $seq_len = 0;
	for (my $i = 1; $i <= $length; $i++) {
		
		if ($seq1->pos($i) ne '-') {
		
			if  ($seq1->pos($i) eq $seq2->pos($i)) {

				$identity++;
			}
			$seq_len++;
		}
	}

	return $identity/$seq_len;
}

###############################################################
# General output methods
###############################################################

#-------------------------------------------------------------#
#
# to_fastq
#
# Description:
#	To convert an object sequence to fastq
#
# Usage:
#	$string = $object->to_fastq();
#
# Arguments:
#	UPPERCASE => 0|1  (default 0)
#   LOWERCASE => 0|1  (default 0)
#   LINELEN   => int  (default 60)
#
# Return:
#	string containing the fasta output
#-------------------------------------------------------------#
sub to_fastq() {
	
	my $self = shift;
	my %opt  = @_;
	my $out;
	my $seq;
    my $qual;
	my $ox = $self->get_ox();
    my $oc = $self->get_oc();
    my $os = $self->get_os();
    my $linelen = (defined $opt{LINELEN})? $opt{LINELEN} : 60;
	$seq  = $self->get_seq();
	$qual = $self->get_qual();
	$seq  = lc($seq) if ($opt{'LOWERCASE'});
	$seq  = uc($seq) if ($opt{'UPPERCASE'});
	$out  = "@".$self->get_id()." ".$self->get_des();
	$out .= " $ox" if ($ox);
	$out .= " TAXONOMY=\"$oc $os\"" if ($oc);
	$out .= "\n";
    if ($linelen > 0)
    {
	    $out .= "$1\n" while ($seq =~ /(\S{1,$linelen})/g); 
	    $out .= "+\n";
        $out .= "$1\n" while ($qual =~ /(\S{1,$linelen})/g);
    }
    else
    {
        $out .= $seq."\n";
	    $out .= "+\n";
        $out .= $qual."\n";
    }
    return $out;
}

#-------------------------------------------------------------#
#
# to_fasta
#
# Description:
#	To convert an object sequence to fasta
#
# Usage:
#	$string = $object->to_fasta();
#
# Arguments:
#	UPPERCASE => 0|1  (default 0)
#   LOWERCASE => 0|1  (default 0)
#   LINELEN   => int  (default 60)
#
# Return:
#	string containing the fasta output
#-------------------------------------------------------------#
sub to_fasta() {
	
	my $self = shift;
	my %opt  = @_;
	my $out;
	my $seq;
	my $ox = $self->get_ox();
    my $oc = $self->get_oc();
    my $os = $self->get_os();
    my $zv = $self->get_zv();
    my $linelen = (defined $opt{LINELEN})? $opt{LINELEN} : 60;
	$seq  = $self->get_seq();
	$seq  = lc($seq) if ($opt{'LOWERCASE'});
	$seq  = uc($seq) if ($opt{'UPPERCASE'});
	$out  = ">".$self->get_id()." ".$self->get_des();
	$out .= " $ox" if ($ox);
	$out .= " TAXONOMY=\"$oc $os\"" if ($oc);
	$out .= " VALID_OS=\"$zv\"" if ($zv);
    $out .= "\n";
    if ($linelen > 0)
    {
	    $out .= "$1\n" while ($seq =~ /(\S{1,$linelen})/g); 
    }
    else
    {
        $out .= $seq."\n";
    }
	return $out;
}

#-------------------------------------------------------------#
#
# to_sprot
#
# Description:
#	To convert an object sequence to swissprot format
#
# Usage:
#	$string = $object->to_sprot();
#
# Return:
#	string containing the swissprot output
#-------------------------------------------------------------#
sub to_sprot() {

	my $self = shift;
	my $out;
	my $seq  = $self->get_seq();
	
	# ID
	my $out = "ID   ".$self->get_id() . "     STANDARD;      ";
	$out .= ($self->get_type() eq 'P')? 'PRT;   ' : 'DNA;   ';
	$out .= $self->get_length();
	$out .= ($self->get_type() eq 'P')? " AA.\n" : " BP.\n";
	
	# AC
	$out .= "AC   ";
	foreach ($self->get_ac()) {
	
		$out   .= "$_ ";
	}
	$out .= "\n";
	
	# DE
	$out   .= "DE   ".$self->get_des()."\n";

	# DR
	foreach ($self->get_dr()) {

		$out .= "DR   $_.\n";
	}
	
    # OC
	foreach ($self->get_oc()) {

		$out .= "OC   $_\n";
	}
    
    # OS
	foreach ($self->get_os()) {

		$out .= "OS   $_\n";
	}

	# FT
	foreach ($self->get_ft()) {

		$out .= "FT   $_.\n";
	}

	# SQ
	$out .= "SQ   SEQUENCE   " . $self->get_length() ." ";
	$out .= ($self->get_type() eq 'P')? 'AA; ' : 'BP; ';
	$out .= $self->get_checksum() . " CRC64;";

	my $i = 0;
	while ($seq =~ /(\S{1,10})/g) {

		if ($i % 6 == 0) {
		
			$out .= "\n     ";
		}
		else {

			$out .= " ";
		}

		$out .= $1;
		$i++;
	}
	chomp($out);
	$out .= "\n//\n";
	return $out;
}

#-------------------------------------------------------------#
#
# to_dat
#
# Description:
#	To convert an object sequence to swissprot format
#
# Usage:
#	$string = $object->to_sprot();
#
# Return:
#	string containing the dat output
#-------------------------------------------------------------#
sub to_dat() 
{
	my $self = shift;
	my $out;
	my $seq  = $self->get_seq();
	
	# ID
	my $out = "ID   ".$self->get_id()."\n";
	
	# AC
	$out .= "AC   ".join(';',$self->get_ac()).";\n";
	
	# DE
	$out   .= "DE   ".$self->get_des()."\n";

	# DR
	foreach ($self->get_dr()) {

		$out .= "DR   $_.\n";
	}
	
    # OC
	foreach ($self->get_oc()) {

		$out .= "OC   $_\n";
	}
    
    # OS
	foreach ($self->get_os()) {

		$out .= "OS   $_\n";
	}

	# FT
	foreach ($self->get_ft()) {

		$out .= "FT   $_.\n";
	}

	# SQ
	$out .= "SQ   SEQUENCE   " . $self->get_length() ." ";
	$out .= ($self->get_type() eq 'P')? 'AA; ' : 'BP; ';
	$out .= $self->get_checksum() . " CRC64;";

	my $i = 0;
	while ($seq =~ /(\S{1,10})/g) {

		if ($i % 6 == 0) {
		
			$out .= "\n     ";
		}
		else {

			$out .= " ";
		}

		$out .= $1;
		$i++;
	}
	chomp($out);
	$out .= "\n//\n";
	return $out;
}

#-------------------------------------------------------------#
#
# to_pir
#
# Description:
#	To convert an object sequence to pir
#
# Usage:
#	$string = $object->to_pir();
#
# Return:
#	string containing the pir output
#-------------------------------------------------------------#
sub to_pir {
	
	my $self = shift;
	my $out;
	my $seq;

	$seq  = $self->get_seq();
	$out  = ">".$self->get_type()."1;".$self->get_id()."\n";
	$out .= $self->get_des()."\n";
	$out .= "$1\n" while ($seq =~ /(\S{1,60})/g); 
	chomp($out);
	$out .= "*\n";
	return $out;
}

#-------------------------------------------------------------#
#
# to_msf
#
# Description:
#	To convert an object sequence to msf
#
# Usage:
#	$string = $object->to_msf();
#
# Return:
#	string containing the msf output
#-------------------------------------------------------------#
sub to_msf {
	
	my $self = shift;
	my $out;
	my $seq;

	$out  = " ".$self->get_id().".msf  ";
	$out .= "MSF: ".$self->get_length()."  ";
	$out .= "Type: ".$self->get_type()."  ";
	$out .= "Check: ".$self->get_checksum()." ";
	$out .= "..\n\n";

	$out .= sprintf(" Name: %15s",$self->get_id());
	$out .= sprintf("\tLen:%4i  ",$self->get_length());
	$out .= sprintf("Check:%6i  ",$self->get_checksum());
	$out .= sprintf("Weight:  %4.2f\n",1);
	$out .= "\n//\n";
	
	$seq  = $self->get_seq();
	
	while ($seq =~ /(\w{1,50})/g) {
		
		my $tmp = $1;

		$out .= "\n";
		$out .= sprintf("%15s  ",$self->get_id());
		$out .= "$1 " while ($tmp =~ /(\S{1,10})/g);
		$out .= "\n";
	}

	return $out;
}

###############################################################
# Method to extract a subsequence or return a midified sequence
###############################################################

#-------------------------------------------------------------#
# subseq()
# 
# Description:
#	To extract a sub sequence 
#
# Usage:
#	$string = $object->subseq(int,int);
#
#	where the two int values represent the start and end positions
#	in the sequence.
#
# Return:
#	string containing the sub-sequence
#-------------------------------------------------------------#
sub subseq() {

	my ($self,$start,$end) = @_;
	my ($offset,$length);
	
	if ($start > $end) {

		my $tmp = $start;
		$start = $end;
		$end = $tmp;
	}

	if ($start < 1 || $end-1 > $self->get_length()) {

		print STDERR "Sequence boundary error in subseq()\n";
		return undef;
	}

	$length = $end-$start+1;
	$offset = $start-1;
		
	return substr($self->get_seq(),$offset,$length);
}

#-------------------------------------------------------------#
#
# reverse()
#
# Description:
#	To reverse a sequence
#
# Usage:
#	$object->reverse();
#
# Return:
# 	self
#-------------------------------------------------------------#
sub reverse {
	my @tmp = split //,$_[0]->get_seq();
	$_[0]->set_seq(join('',reverse(@tmp)));
    @tmp = ();
    @tmp = split //,$_[0]->get_qual();
    $_[0]->set_qual(join('',reverse(@tmp)));
    return $_[0];
}

#-------------------------------------------------------------#
#
# complement()
#
# Description:
#	To complement a DNA sequence
#
# Usage:
#	$object->complement();
#
# Return:
# 	self
#-------------------------------------------------------------#
sub complement {

	my $i;
	my @tmp = split //,$_[0]->get_seq();
	
	for($i=0;$i<scalar(@tmp);$i++) {

		$tmp[$i] = $SeqCodes::complement{$tmp[$i]};
	}
	
	$_[0]->set_seq(join('',@tmp));
    return $_[0];
}


#-------------------------------------------------------------#
#
# window_reverse()
#
# Description:
#	To reverse a sequence in windows of lengrh x
#
# Usage:
#	$object->window_reverse(x);
#
# Arguments:
#	int x: size of the window to permute
#
# Return:
#	void
#
#-------------------------------------------------------------#
sub window_reverse() {

	my ($self,$window) = @_;
	my ($i,@seq,$seq,$final_seq);
	
	$seq = $self->get_seq();
	push(@seq,$1) while ($seq =~ /(.{1,$window})/g);
	
	for ($i=0;$i<scalar(@seq);$i+=2) {
			
		$final_seq .= $seq[$i+1].$seq[$i];
	}

	$self->set_seq($final_seq);
    return $self;
}

#-------------------------------------------------------------#
#
# window_shuffle()
#
# Description:
#	To shuffle a sequence in windows of lengrh x
#
# Usage:
#	$object->window_shuffle(x);
#
# Arguments:
#	int x: size of the window to shuffle
#
# Return:
#	void
#
#-------------------------------------------------------------#
sub window_shuffle() {

	my ($self,$window) = @_;
	my @seq;
	my $seq = $self->get_seq();

	push(@seq,$1) while ($seq =~ /(.{1,$window})/g);
	
	for (my $i = 0; $i < scalar(@seq); $i++) {
		
		my @tmp  = split //,$seq[$i];
        array_fisher_yates_shuffle(\@tmp);
		$seq[$i] = join('',@tmp);
	}
	
	$self->set_seq(join('',@seq));
	return $self;
}

#-------------------------------------------------------------#
#
# permute()
#
# Description
#	To permute a sequence based on windows of size x
#
# Usage:
#	$object->permute()
#
# Arguments:
#	int x: size of the windows
#
# Return:
#	void
#
#-------------------------------------------------------------#	
sub permute() {
	
	my ($self,$window) = @_;
	my ($i,@seq,$seq,$final_seq);

	$seq = $self->get_seq();
	push(@seq,$1) while ($seq =~ /(.{1,$window})/g);
    array_fisher_yates_shuffle(\@seq);
	$final_seq = join('',@seq);
	
	$self->set_seq($final_seq);
    return $self;
}


###############################################################
# Methods to modify sequences
###############################################################

#-------------------------------------------------------------#
#
# set_seq()
# 	To set a sequence
#
# Usage:
#	$object->set_seq($string);
# 
# Arguments:
#	$string: The sequnce
# 
# Return
#	Void
#
#-------------------------------------------------------------#
sub set_seq() {

	my ($self,$seq) = @_;

	$self->set_length(length($seq));
	$self->{'seq'} = $seq;
}

###############################################################
# FASTQ
package fastq;

@fastq::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
#
# read()
#
# Description:
#	Method to read a fastq entry into an object. Uses the to_object
#	method.
#
# Usage:
#	$fastq->read(\*FILENAME)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub read() {

	my ($self,$fh) = @_;
	
	return $self->to_object($fh);
}

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to read a fasta entry into an object
#
# Usage:
#	$fasta->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object() {

	my ($self,$fh) = @_;
	my $offset = 0;
	my $sequence;
    my $quality;
	my $flag;
	
	if (!defined($fh)) {
	
		confess("ERROR:",ref($self),"->to_object() called without filehandle\n");
	}

	$self->clear();
	
	LOOP: while (!eof($fh)) {
		my $c = BioSeqFormat::peek($fh); 
		chomp($c);
        if ($c eq '@' && $flag!=2) {
			
			if (defined($self->get_id())) {
				
				seek($fh,$offset,0) if ($offset > -1); # Only if we are not reading from a pipeline!;
				last LOOP;
			}
			
			<$fh> =~ /^\@\s*(\S+)(?:\s+(.+))?\s*$/;
			$self->set_id($1);
			$self->set_des($2) if ($2);
			$flag = 1;
		}
		if ($c eq '+' && 1 == $flag) 
        {
            $flag = 2;
            <$fh>;
		}
		elsif ( $c eq ' ' || $c eq '') {

			<$fh>; # Skip empty lines and html tags
		}
		elsif ($flag) {
			
			my $tmp = <$fh>;
			chomp($tmp);
            if (1 == $flag)
            {   
                $sequence .= $tmp;
		    }
            elsif (2 == $flag) 
            {
                $quality .= $tmp;
                $flag = 1 if (length($sequence) == length($quality));
            }
        }
		
		$offset = tell;
	}

	$self->set_type('N');
	$self->set_checksum(BioSeqFormat::checksum($sequence));
	$self->set_seq($sequence);
    $self->set_qual($quality);
    return ($flag)? 1:undef;
}



###############################################################
# FASTA
package fasta;

@fasta::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
#
# read()
#
# Description:
#	Method to read a fasta entry into an object. Uses the to_object
#	method.
#
# Usage:
#	$fasta->read(\*FILENAME)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub read() {

	my ($self,$fh) = @_;
	
	return $self->to_object($fh);
}

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to read a fasta entry into an object
#
# Usage:
#	$fasta->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object() {

	my ($self,$fh) = @_;
	my $offset = 0;
	my $sequence;
	my $flag;
	
	if (!defined($fh)) {
	
		confess("ERROR:",ref($self),"->to_object() called without filehandle\n");
	}

	$self->clear();
	
	LOOP: while (!eof($fh)) {
		
		my $c = BioSeqFormat::peek($fh); 
		chomp($c);
		
		if ($c eq '>') {
			
			if (defined($self->get_id())) {
				
				seek($fh,$offset,0) if ($offset > -1); # Only if we are not reading from a pipeline!;
				last LOOP;
			}
			
			<$fh> =~ /^>\s*(\S+)(?:\s+(.+))?\s*$/;
			$self->set_id($1);
			$self->set_des($2) if ($2);
			$flag = 1;
		}
		elsif ($c eq '<' || $c eq ' ' || $c eq '') {

			<$fh>; # Skip empty lines and html tags
		}
		elsif ($flag) {
			
			my $tmp    = <$fh>;
			$tmp       =~ s/\s+//g;
			$sequence .= $tmp;
		}
		
		$offset = tell;
	}

	$self->set_type(BioSeqFormat::guess_type($sequence));
	$self->set_checksum(BioSeqFormat::checksum($sequence));
	$self->set_seq($sequence);
	return ($flag)? 1:undef;
}

###############################################################
# DAT (ENA)
package dat;

@dat::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to read a fasta entry into an object
#
# Usage:
#	$dat->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object {

	my ($self,$fh) = @_;
	my $flag;

	$self->clear();	

	while (<$fh>) {
		
		if ($flag) {

			if (/^\/\//) 
            {
				last;
			}
			
			if (/^AC\s+(\S+)/) {
				
				$self->add_ac(split /;/,$1);
				next;
			}
			
            if (/^DE\s+(\S.*)$/) {

				if ($self->get_des()) 
                {
					$self->set_des($self->get_des()." ".$1);
				}
				else 
                {
					$self->set_des($1);
				}
				next;
			}
			
			if (/^OS\s+(\S.*)$/) {

				if ($self->get_os()) 
                {
					$self->set_os($self->get_os()." ".$1);
				}
				else 
                {
					$self->set_os($1);
				}
				next;
			}

			if (/^OC\s+(\S.*)$/) 
            {
				if ($self->get_oc()) 
                {	
					$self->set_oc($self->get_oc()." ".$1);
				}
				else 
                {	
                    $self->set_oc($1);
				}
				next;
			}
			
			if (/^OX\s+(\S.*)$/) 
            {
				$self->set_ox($1);
				next;
			}

			if (/^DR\s+(\S.+)\.$/) 
            {	
				$self->add_dr($1);
				next;
			}
			
			if (/^FT\s+(\S.+)\.$/) 
            {	
				$self->add_ft($1);
				next;
			}
			
            if (/^\s+([\D]+)\s*(\d+)?\s*$/)
			{
				my $tmp = $1;
				$tmp =~ s/\s+//g;
				$self->set_seq($self->get_seq().$tmp);
				next;
			}
		}		
		
		if (/^ID\s+(\S.*)/) 
        {
			$self->set_id($1);
			$flag = 1;
		}
	}
	
	$self->set_type(BioSeqFormat::guess_type($self->get_seq()));

	$self->set_checksum(BioSeqFormat::checksum($self->get_seq()));
	
    return ($self->get_id())? 1:undef;
}


###############################################################
# SPROT_TREMBL
package sprot;

@sprot::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to read a fasta entry into an object
#
# Usage:
#	$sprot->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object {

	my ($self,$fh) = @_;
	my $flag;

	$self->clear();	

	while (<$fh>) {
		
		if ($flag) {

			if (/^\/\//) {

				last;
			}
			
			if (/^AC\s+(\S+)/) {
				
				$self->add_ac(split /;/,$1);
				next;
			}

			if (/^DE\s+(\S.*)$/) {

				if ($self->get_des()) {
				
					$self->set_des($self->get_des()." ".$1);
				}
				else {
					$self->set_des($1);
				}
				next;
			}
			
			if (/^OS\s+(\S.*)$/) {

				if ($self->get_os()) {
				
					$self->set_os($self->get_os()." ".$1);
				}
				else {
					$self->set_os($1);
				}
				next;
			}

			if (/^OC\s+(\S.*)$/) {

				if ($self->get_oc()) {
				
					$self->set_oc($self->get_oc()." ".$1);
				}
				else {
					$self->set_oc($1);
				}
				next;
			}
			
			if (/^OX\s+(\S.*)$/) {

				$self->set_ox($1);
				next;
			}

			if (/^DR\s+(\S.+)\.$/) {
				
				$self->add_dr($1);
				next;
			}
			
			if (/^FT\s+(\S.+)\.$/) {
				
				$self->add_ft($1);
				next;
			}
			#if (/^\s+(\D+)\s*(\d+)?\s*$/ || /^SQ\s+(\S+)/) {
			if (/^\s+(\D+)\s*(\d+)?\s*$/)
			{
				my $tmp = $1;
				$tmp =~ s/\s+//g;
				$self->set_seq($self->get_seq().$tmp);
				next;
			}
		}		
		
		if (/^ID\s+(\S+)/) {

			$self->set_id($1);
			$flag = 1;
		}
		
	}
	
	$self->set_type(BioSeqFormat::guess_type($self->get_seq()));
	$self->set_checksum(BioSeqFormat::checksum($self->get_seq()));
	return ($self->get_id())? 1:undef;
}


###############################################################
# PIR
package pir;

@pir::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to read a pir format into an object
#
# Usage:
#	$pir->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object {

	my ($self,$fh) = @_;
	my $offset = 0;

	$self->clear();
	
	while (<$fh>) {
		
		# Read sequence
		if (/^\s*([^\:\;\>]+)/) {

			my $tmp = $1;
		
			$tmp =~ s/\s+//g;
			$tmp =~ s/\*$//;
			$self->set_seq($self->get_seq().$tmp);
		}
		
		# Read '>' line
		elsif (/^>[PF\d]+;\s*(\S+)/) {

			if ($self->get_id()) {
				
				seek($fh,$offset,0);
				last;
			}
			
			$self->clear();
			$self->set_id($1);
			$_ = <$fh>; 
			chomp;
			$self->set_des($_);
		}
		
		$offset = tell;
	}

	$self->set_type(BioSeqFormat::guess_type($self->get_seq()));
	return ($self->get_id())? 1:undef;
}

###############################################################
# QDAT
package qdat;

@qdat::ISA = qw(BioSeqFormat);

#-------------------------------------------------------------#
# 
# to_object
#
# Description:
#	Method to qdat a fasta entry into an object
#
# Usage:
#	$qdat->to_object(\*FILEHANDLE)
#
# Return:
#	1 if ok, otherwise undef
#
#-------------------------------------------------------------#
sub to_object {

    my ($self,$fh) = @_;
    my $flag;

    $self->clear();	

    while (<$fh>) {

        if ($flag) {

            if (/^\/\//) {

                last;
            }

            if (/^AC\s+(\S+)/) {

                $self->add_ac(split /;/,$1);
                next;
            }

            if (/^DE\s+(\S.*)$/) {

                if ($self->get_des()) {

                    $self->set_des($self->get_des()." ".$1);
                }
                else {
                    $self->set_des($1);
                }
                next;
            }

            if (/^OS\s+(\S.*)$/) {

                if ($self->get_os()) {

                    $self->set_os($self->get_os()." ".$1);
                }
                else {
                    $self->set_os($1);
                }
                next;
            }

            if (/^OC\s+(\S.*)$/) {

                if ($self->get_oc()) {

                    $self->set_oc($self->get_oc()." ".$1);
                }
                else {
                    $self->set_oc($1);
                }
                next;
            }
            
            if (/^ZV\s+(\S.*)$/) 
            {
                $self->set_zv($1);
                next;
            }

            if (/^OX\s+(\S.*)$/) {

                $self->set_ox($1);
                next;
            }

            if (/^DR\s+(\S.+)\.$/) {

                $self->add_dr($1);
                next;
            }

            if (/^FT\s+(\S.+)\.$/) {

                $self->add_ft($1);
                next;
            }
            if (/^SQ\s+(\S+)/) {
                {
                    my $tmp = $1;
                    $tmp =~ s/\s+//g;
                    $self->set_seq($self->get_seq().$tmp);
                    next;
                }
            }		
        }

        if (/^ID\s+([^;\s]+)/) {

            $self->set_id($1);
            $flag = 1;
        }

    }

    $self->set_type(BioSeqFormat::guess_type($self->get_seq()));
    $self->set_checksum(BioSeqFormat::checksum($self->get_seq()));
    return ($self->get_id())? 1:undef;
}

1;


