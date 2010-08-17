# $Id: Primer3.pm,v 0.01 2007-03-27 12:43:27 heikki Exp $
#
# BioPerl module for Bio::Tools::Primer3::Sequence
#
# Cared for by Chris Fields cjfields at bioperl dot org
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer3::Sequence - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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

package Bio::Tools::Primer3::Sequence;

use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new
 Function: Builds a new Bio::Tools::Primer3::Sequence object 
 Returns : an instance of Bio::Tools::Primer3::Sequence
 Args    : 
 
=cut

sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
	my ($seq) = $self->_rearrange();
	return $self;
}

1;
