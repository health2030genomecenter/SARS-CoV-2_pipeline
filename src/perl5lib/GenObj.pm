package GenObj;
$VERSION = 0.4.1;

#-------------------------------------------------------------#	
###############################################################	
###############################################################
#
# This module is developed and maintained by Lorenzo CERUTTI
# (lorenzo.cerutti@lulix.net)
#
# You can use it and modify it, just leave this comment and 
# author information. ... be nice with the author :)
#
###############################################################
###############################################################
#-------------------------------------------------------------#	


###############################################################
# DEBUG HISTORY
#
# VERSION 0.4.1
#	Add a general constant to set tmp directory.
#
# VERSION 0.4
#	Now uses Carp module and confess if an error is detect.
#	'Set' methods now can be used for list attributes passing:
#		first: the position in the list
#		second: the value
#		ex: $self->set_xxx(2,$value);
#
# VERSION 0.3
#	clone(): Solved cloning of hash
#	         "Can't call method "clone" on unblessed reference at
#            GenObj.pm line 443"
#
###############################################################	

use constant DEBUG => undef;
#use constant DEBUG => 1;

use constant LIST	=> 'list';
use constant VAR	=> 'var';
use constant TMP    => '/tmp/';
use Carp;

use strict;
no strict "refs";
use vars '$AUTOLOAD';

###############################################################
#
# GenaralObject.pm
#
#   Package to generate dinamically general methods for objects.
#
# Lorenzo Cerutti, 10/4/2001
#
###############################################################

###############################################################
# 
# new() 
#	General constructor for an object.
#
# Usage:
#	$object = CLASS->new();
#
# Arguments:
#	CLASS: class inheriting GeneralObject
#
# Return:
#	$object: instance of CLASS
#
###############################################################
sub new { 
	
	my $self = bless {}, ref($_[0]) || $_[0];
	
	foreach ($self->get_attributes()) {
		
		my ($attr_name,$attr_type) = split /\:/;

		print STDERR "GenObj::new() ->($attr_name,$attr_type)<\n" if (DEBUG);
		
		CASE: {
			
			if ($attr_type eq 'val')  {
				$self->{$attr_name} = undef; 
				last CASE; 
			}
			if ($attr_type eq 'list') { 
				$self->{$attr_name} = ();
				last CASE; 
			}
			confess("ERROR: \"$attr_name\" has an error in its type \"$attr_type\"");
		}
	}

	return $self;
}

###############################################################
#
# is_defined()
#	Return true if the attribute has already been declared
#
# Usage:
#	$bool = $object->is_defined($attribute);
#
# Arguments:
#	$attribute: name of the attribute to check
#
# Return:
#	1 if attribute is declared, else UNDEF
#
################################################################
sub is_defined {

	return ($_[0]{$_[1]})? 1:undef;
}

###############################################################
#
# get_parent_names()
#	Return an array of direct parent classes
#
# Usage:
#	@list = $object->get_parent_names()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing the list of direct parents of the 
#			object
#
###############################################################
sub get_parent_names {

	my $self = ref($_[0]) || $_[0];
	
	return @{"${self}::ISA"};
}

###############################################################
#
# get_all_parent_names()
#	Return an array of parent classes recursivelly
#
# Usage:
#	@list = $object->get_all_parent_names()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing the list of parents of the object
#
###############################################################
sub get_all_parent_names {

	my $self = ref($_[0]) || $_[0];
	
	my @parents;
	
	if (@{"${self}::ISA"}) {

		foreach my $parent_class ($self->get_parent_names()) {

			if (defined($parent_class->get_parent_names())) {

				push @parents, get_all_parent_names($parent_class);
			}
			
			push @parents, $parent_class;
		}
	}
	return @parents;
}

###############################################################
#
# is_parent()
#	Tell if a name is the direct parent of a class
#
# Usage:
#	$bool = $object->is_parent($string)
#
# Arguments:
#	$string    Name of the parent
#
# Return:
#	1 if the $string correspond to the parent name
#	undef otherwise
#
###############################################################
sub is_parent {

	my ($self,$name) = @_;

	foreach my $parent_name ($self->get_parent_names()) {

		return 1 if ($name eq $parent_name);
	}
	
	return undef;
}

###############################################################
#
# get_attributes()
#	Return an array of attributes of the object
#
# Usage:
#	@list = $object->get_attribute_names()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing attributes of the object
#
###############################################################
sub get_attributes {

	my $self = ref($_[0]) || $_[0];
	
	return @{"${self}::_ATTRIBUTES_"};
}

###############################################################
#
# get_attribute_names()
#	Return an array of attribute names of the object
#
# Usage:
#	@list = $object->get_attribute_names()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing attribute names of the object
#
###############################################################
sub get_attribute_names {

	my $self = ref($_[0]) || $_[0];
	my @attr_names;

	foreach ($self->get_attributes()) {
	
		my ($name,$type) = split /\:/;
		push(@attr_names,$name);
	}
	
	return @attr_names;
}

###############################################################
#
# get_all_attributes()
#	Return an array of attributes redundantly
#
# Usage:
#	@list = $object->get_all_attributes()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing attributes of the object and
#		its parents
#
###############################################################
sub get_all_attributes {

	my $self = ref($_[0]) || $_[0];

	my @attributes = $self->get_attributes();
	
	if (defined($self->get_parent_names())) {

		foreach my $parent_class ($self->get_parent_names()) {

			push @attributes, get_all_attributes($parent_class);
		}
	}
	return @attributes;
}


###############################################################
#
# get_all_attribute_names()
#	Return an array of attributes redundantly
#
# Usage:
#	@list = $object->get_all_attribute_names()
#
# Arguments:
#	void
#
# Return:
#	@list	array containing attributes of the object and is
#			its parents
#
###############################################################
sub get_all_attribute_names {

	my $self = ref($_[0]) || $_[0];
	
	print STDERR "CALL: GenObj::get_all_attribute_names() from $self\n" if DEBUG;
	
	my @attributes = $self->get_attribute_names();
	
	if (defined($self->get_parent_names())) {

		foreach my $parent_class ($self->get_parent_names()) {

			push @attributes, get_all_attribute_names($parent_class);
		}
	}
	return @attributes;
}

###############################################################
#
# which_type()
#	Return the type of an attribute
#
# Usage:
#	$type = $object->which_type($string)
#
# Arguments:
#	$string: Name of the attribute to check
#
# Return:
#	$string: Type of the attribute
#
###############################################################
sub which_type {

	my ($self,$name) = @_;

	foreach ($self->get_all_attributes()) {

		my ($attr,$type) = split /\:/;
		
		return $type if ($attr eq $name);
	}

	return undef;
}

###############################################################
#
# is_attribute()
#	Check if a name is a direct attribute of the class
#
# Usage:
#	$bool = $object->is_attribute($string)
#
# Arguments:
#	$string		Name on the attribute to check
#
# Return:
#	1 if $string is an attribute
#	undef otherwise
#
###############################################################
sub is_attribute {

	my ($self,$name) = @_;

	foreach my $attribute ($self->get_attribute_names()) {

		return 1 if ($attribute eq $name);
	}

	return undef;
}

###############################################################
#
# add_attribute_name()
#	Add an attribute name to the class
#
# Usage:
#	$object->add_attribute_name(@names)
#
# Arguments:
#	@names	list of attributes names
#
# Return:
#	void
#
###############################################################
sub add_attribute_name {
	
	my ($self, @names) = @_;
	
	$self = ref($self) || $self;
	
	push @{"${self}::_ATTRIBUTES_"}, @names;
}


###########################################################
#
# clear()
#	General method to clear an object.
#
# Usage:
#	$object->clear()
#
# Arguments:
#	void
#
# Return:
#	void
#
###############################################################
sub clear {
	
	my $self = shift;
	
	foreach my $attribute ($self->get_all_attribute_names()) {
		
		my $class = ref($self->{$attribute});
		$class = ($class eq "ARRAY" || $class eq "HASH")? undef:$class;
		
		$self->{$attribute}->clear() if (ref($class));

		delete $self->{$attribute}; 
	} 
}

###############################################################
#
# clone()
#	General method to clone an object
#
# Usage:
#	$newObject = $object->clone();
#
# Arguments:
#	void
#
# Return:
#	$object
#
###############################################################
sub clone {

	my $self = shift;
	my $new = bless {}, ref($self) || $self;
	
	foreach my $attribute ($self->get_all_attribute_names()) {
		
		if ( ref($self->{$attribute} ) eq "ARRAY") {

			foreach my $value (@{$self->{$attribute}}) {
		
				if ( ref($value) ) {
					
					push @{$new->{$attribute}}, $value->clone();

				}
				else {
					
					push @{$new->{$attribute}}, $value;
				}
			}
		}
		elsif ( ref($self->{$attribute}) eq 'HASH' ) {
			
			$new->{$attribute} = deepcopy($self->{$attribute});
		}

		elsif ( ref($self->{$attribute}) ) {
			
			$new->{$attribute} = $self->{$attribute}->clone();
		}
		else {
			
			$new->{$attribute} = $self->{$attribute};
		}
	}
	return $new;
}

sub deepcopy {

	if (ref $_[0] eq 'HASH') {
		
		return { map(deepcopy($_), %{$_[0]}) };
	} 
	elsif (ref $_[0] eq 'ARRAY') {

		return [ map(deepcopy($_), @{$_[0]}) ];
	}
	$_[0];
}

###############################################################
#
# AUTOLOAD METHODS
#
# set_..(), add_..(), get_..(), getstring_..(),  getref_..(), 
# push_..(), shift_..(), delete_..()
#
# Define a number of general methods.
#
###############################################################
sub AUTOLOAD {
	
	my ($self, @value) = @_;
	
	#########################################################
	# Set methods
	if ($AUTOLOAD =~ /\:set_(\w+)/) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {
				
				my $curr_type = $self->which_type($attribute);
				
				if ($curr_type eq 'val') {

					*{$AUTOLOAD} = sub { 
					
						return $_[0]->{$attribute} = $_[1]; 
					};
		
					return $self->{$attribute} = $value[0];
				}
				else {

					*{$AUTOLOAD} = sub { 
					
						return $_[0]->{$attribute}->[$_[1]] = $_[2]; 
					};
		
					return $self->{$attribute}->[$value[0]] = $value[1];
				}
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}	
	
	#######################################################
	# Push methods
	if ($AUTOLOAD =~ /\:push_(\w+)/) {
		
		my $attribute = $1;
	
		foreach ($self->get_all_attribute_names()) {
		
			if ($attribute eq $_) {
	
				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': push_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { 
					
					my $self = shift;

					return push(@{$self->{$attribute}},@_);
				};
		
				return push(@{$self->{$attribute}},@value); 
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}

	#########################################################
	# Pop methods
	if ( $AUTOLOAD =~ /\:pop_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {

				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': shift_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { return pop(@{$_[0]->{$attribute}}) };
		
				return pop(@{$self->{$attribute}});
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}

	#######################################################
	# Add methods
	if ($AUTOLOAD =~ /\:add_(\w+)/) {
		
		my $attribute = $1;
	
		foreach ($self->get_all_attribute_names()) {
		
			if ($attribute eq $_) {
				
				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': add_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { 
					
					my $self = shift;

					return push(@{$self->{$attribute}},@_);
				};
		
				return push(@{$self->{$attribute}},@value); 
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}
	

	#########################################################
	# Get methods
	if ( $AUTOLOAD =~ /\:get_(\w+)/ ) {
		
		my $attribute = $1;
	
		print STDERR "GenObj::get_$attribute ->$self<-\n" if (DEBUG);

		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {
				
				if ($self->which_type($attribute) eq 'list') { 
			
					
					*{$AUTOLOAD} = sub { 
						
						if (defined($_[1])) {
							
							return ${$_[0]->{$attribute}}[$_[1]];
						}
						else {
							return @{$_[0]->{$attribute}}; 
						}
					};
					
					
					if (defined($value[0])) {
						
						return ${$self->{$attribute}}[$value[0]];

					}
					else {
						return @{$self->{$attribute}};
					}
				}
				else {
		
					*{$AUTOLOAD} = sub { return $_[0]->{$attribute} };
			
					return $self->{$attribute};
				}
			}
		}
			
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}
	
	#########################################################
	# Get2string methods
	if ( $AUTOLOAD =~ /\:getstring_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {
				
				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': getstring_$attribute not possible: attribute must be of type 'list':");
				}

				my $separator = ($value[0])?$value[0]:'';
				
				*{$AUTOLOAD} = sub { 
					
					my $s = ($_[1])?$_[1]:'';
					return join($s,@{$_[0]->{$attribute}}) 
				};
		
				return join($separator,@{$self->{$attribute}});
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}
	
	#########################################################
	# Getref methods
	if ( $AUTOLOAD =~ /\:getref_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {

				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': getref_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { 
				
					return \@{$_[0]->{$attribute}};
				};
		
				return \@{$self->{$attribute}};
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}
	
	#########################################################
	# Shift methods
	if ( $AUTOLOAD =~ /\:shift_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {

				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': shift_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { return shift(@{$_[0]->{$attribute}}) };
		
				return shift(@{$self->{$attribute}});
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}

	#########################################################
	# Unshift methods
	if ( $AUTOLOAD =~ /\:unshift_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {

				if ($self->which_type($attribute) ne 'list') {
				
					confess("ERROR executing '$AUTOLOAD': unshift_$attribute not possible: attribute must be of type 'list':");
				}

				*{$AUTOLOAD} = sub { 
					
					my ($self,@list) = @_;
					return unshift(@{$_[0]->{$attribute}},@list);
				};
		
				return unshift(@{$self->{$attribute}},@value);
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}

	#########################################################
	# Delete methods
	if ( $AUTOLOAD =~ /\:delete_(\w+)/ ) {
		
		my $attribute = $1;
		
		foreach ($self->get_all_attribute_names()) {

			if ($attribute eq $_) {
			
				if ($self->which_type($attribute) eq 'list') {
				
					*{$AUTOLOAD} = sub { return @{$_[0]->{$attribute}} = () };
			
					return @{$self->{$attribute}} = ();
				}
				else {
			
					*{$AUTOLOAD} = sub { return $_[0]->{$attribute} = undef };
			
					return $self->{$attribute} = undef;
				}
			}
		}
		
		confess("ERROR executing '$AUTOLOAD': attribute $attribute not defined for class ",ref($self));
	}
}

1;
