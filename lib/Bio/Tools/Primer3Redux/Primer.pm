# $Id: Report.pm,v 0.01 2007-03-27 12:43:27 heikki Exp $
#
# BioPerl module for Bio::Tools::Primer3Redux::Primer
#
# Cared for by Chris Fields cjfields at bioperl dot org
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer3Redux::Primer - Simple Decorator of a Bio::SeqFeature::Generic
with convenience methods for retrieving Tm, GC, validating primer seq against
attached sequence, etc.

=head1 SYNOPSIS

 # get the Bio::Tools::Primer3Redux::Primer through Bio::Tools::Primer3Redux...
 
 # dies with an error if no sequence is attached, or if sequence region
 # does not match cached sequence from Primer3.  Useful if decorating an already
 # generated Bio::Seq with primers.
 
 $primer->validate_seq;
 
 my $seq = $primer->seq; # Bio::Seq object
 if ($primer->melting_temp < 55) {
	warn "Primer ".$primer->display_name." is below optimal temp";
 }
 
 # if primer3 EXPLAIN settings are used...
 print "Run parameters:".$primer->run_description."\n";
 
=head1 DESCRIPTION



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

package Bio::Tools::Primer3Redux::Primer;

use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SeqFeature::Generic);

=head2 oligo_type

 Title    : oligo_type
 Usage    : $obj->oligo_type
 Function : 
 Returns  : 
 Args     : 

=cut

sub oligo_type {
	my ($self, $type) = @_;
	if (defined $type) {
		$self->remove_tag('type') if $self->has_tag('type');
		$self->add_tag_value('type', $type);
	}
	$self->has_tag('type') ? return ($self->get_tag_values('type'))[0] : return;
}

=head2 validate_seq

 Title    : validate_seq
 Usage    : $obj->validate_seq
 Function : 
 Returns  : 
 Args     : 

=cut

sub validate_seq {
	my ($self) = shift;
	my $cached = $self->has_tag('sequence') ? ($self->get_tag_values('sequence'))[0] : '';
	my $seq = $self->seq->seq;
	if ($cached ne $seq) {
		$self->warn("Sequence [$seq] does not match predicted [$cached], check attached sequence");
		return 0;
	}
	return 1;
}

=head2 melting_temp

 Title    : melting_temp
 Usage    : $obj->melting_temp
 Function : 
 Returns  : 
 Args     : 

=cut

sub melting_temp {
	my ($self, $tm) = @_;
	if (defined $tm) {
		$self->remove_tag('tm') if $self->has_tag('tm');
		$self->add_tag_value('tm', $tm);
	}
	$self->has_tag('tm') ? return ($self->get_tag_values('tm'))[0] : return;
}

=head2 gc_content

 Title    : gc
 Usage    : $obj->gc
 Function : 
 Returns  : 
 Args     : 

=cut

sub gc_content {
	my ($self, $gc) = @_;
	if (defined $gc) {
		$self->remove_tag('gc_percent') if $self->has_tag('gc_percent');
		$self->add_tag_value('gc_percent', $gc);
	}
	$self->has_tag('gc_percent') ? return ($self->get_tag_values('gc_percent'))[0] : return;
}

=head2 run_description

 Title    : run_description
 Usage    : $obj->run_description
 Function : 
 Returns  : 
 Args     : 

=cut

sub run_description {
	my ($self, $desc) = @_;
	if (defined $desc) {
		$self->remove_tag('explain') if $self->has_tag('explain');
		$self->add_tag_value('explain', $desc);
	}
	$self->has_tag('explain') ? return ($self->get_tag_values('explain'))[0] : return;
}

1;
