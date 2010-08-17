# $Id: PrimerPair.pm,v 0.01 2007-03-27 12:43:27 heikki Exp $
#
# BioPerl module for Bio::Tools::Primer3Redux::PrimerPair
#
# Cared for by Chris Fields cjfields at bioperl dot org
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer3Redux::PrimerPair - Simple Decorator of a
Bio::SeqFeature::Generic with convenience methods for retrieving left and
right primers, internal oligos, and any amplicon-related information

=head1 SYNOPSIS



=head1 DESCRIPTION

Bio::Tools::Primer3Redux::PrimerPair acts as a simple SeqFeature that bundles
primer pair data together into one object.  This object can be used to retrieve
the amplicon sequence, the forward/reversion (left/right) primers, and any
internal oligos.  Furthermore, any primer information relative to the product
is included as SeqFeature tags.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists
  
=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/
  
=head1 AUTHOR - Chris Fields

  Email cjfields at bioperl dot org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Primer3Redux::PrimerPair;

use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SeqFeature::Generic);

=head2 left_primer

 Title    : left_primer
 Usage    : $obj->left_primer
 Function : 
 Returns  : 
 Args     : 

=cut

sub left_primer {
	shift->forward_primer(@_);
}

=head2 forward_primer

 Title    : forward_primer
 Usage    : $obj->forward_primer
 Function : 
 Returns  : 
 Args     : 

=cut

sub forward_primer {
	my ($self, $primer) = @_;
	if ($primer) {
		$self->throw("Not a Primer object") unless $primer->isa('Bio:::Tools::Primer3Redux::Primer');
		$self->add_SeqFeature($primer, 'EXPAND');
	}
	my ($for) = grep {$_->primary_tag eq 'forward_primer'} $self->get_SeqFeatures;
	return $for;
}

=head2 right_primer

 Title    : right_primer
 Usage    : $obj->right_primer
 Function : 
 Returns  : 
 Args     : 

=cut

sub right_primer { shift->reverse_primer(@_)}

=head2 reverse_primer

 Title    : reverse_primer
 Usage    : $obj->reverse_primer
 Function : 
 Returns  : 
 Args     : 

=cut

sub reverse_primer {
	my ($self, $primer) = @_;
	if ($primer) {
		$self->throw("Not a Primer object") unless $primer->isa('Bio:::Tools::Primer3Redux::Primer');
		$self->add_SeqFeature($primer, 'EXPAND');
	}
	my ($rev) = grep {$_->primary_tag eq 'reverse_primer'} $self->get_SeqFeatures;
	return $rev;
}

=head2 internal_oligo

 Title    : internal_oligo
 Usage    : $obj->internal_oligo
 Function : 
 Returns  : 
 Args     : 

=cut

sub internal_oligo {
	my ($self, $primer) = @_;
	if ($primer) {
		$self->throw("Not a Primer object") unless $primer->isa('Bio:::Tools::Primer3Redux::Primer');
		# Note this doesn't expand to fit; the assumption is this is added
		# after forward/reverse primers are added and acts to ensure the
		# oligo is actually internal to the fragment (otherwise it throws)
		$self->add_SeqFeature($primer);
	}
	my ($oligo) = grep {$_->primary_tag eq 'ss_oligo'} $self->get_SeqFeatures;
	return $oligo;	
}

1;
