# ABSTRACT: Simple Decorator of a Bio::SeqFeature::Generic with convenience
# methods for retrieving Tm, GC, validating primer seq against attached
# sequence, etc.
# AUTHOR:   Chris Fields <cjfields@cpan.org>
# OWNER:    2006-2016 Chris Fields
# LICENSE:  Perl_5

package Bio::Tools::Primer3Redux::Primer;

use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SeqFeature::Generic);
our $AUTOLOAD;

sub oligo_type {
    shift->primary_tag(@_);
}

sub rank {
    my ($self, $rank) = @_;
    if (defined $rank) {
        $self->remove_tag('rank') if $self->has_tag('rank');
        $self->add_tag_value('rank', $rank);
    }
    $self->has_tag('rank') ? return ($self->get_tag_values('rank'))[0] : return;
}

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

sub melting_temp {
    my ($self, $tm) = @_;
    if (defined $tm) {
        $self->remove_tag('tm') if $self->has_tag('tm');
        $self->add_tag_value('tm', $tm);
    }
    $self->has_tag('tm') ? return ($self->get_tag_values('tm'))[0] : return;
}

sub gc_content {
    my ($self, $gc) = @_;
    if (defined $gc) {
        $self->remove_tag('gc_percent') if $self->has_tag('gc_percent');
        $self->add_tag_value('gc_percent', $gc);
    }
    $self->has_tag('gc_percent') ? return ($self->get_tag_values('gc_percent'))[0] : return;
}

sub run_description {
    my ($self, $desc) = @_;
    if (defined $desc) {
        $self->remove_tag('explain') if $self->has_tag('explain');
        $self->add_tag_value('explain', $desc);
    }
    $self->has_tag('explain') ? return ($self->get_tag_values('explain'))[0] : return;
}

sub AUTOLOAD {
  my $self = shift;
  my $type = ref($self) or $self->throw("$self is not an object");
  my $name = $AUTOLOAD;
  $name =~ s/.*://;
  $self->has_tag($name) ? return ($self->get_tag_values($name))[0] : return;
} # AUTOLOAD

# useful for older perl versions
sub DESTROY {}

1;

__END__

# BioPerl module for Bio::Tools::Primer3Redux::Primer
#
# Cared for by Chris Fields cjfields at cpan dot org
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

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

This class is a simple subclass of Bio::SeqFeature::Generic that adds
convenience accessor methods for primer-specific data, such as Tm, GC content,
and other interesting bits of information returned from Primer3.  Beyond that,
the data can be persisted just as any Bio::SeqFeatureI; it doesn't add any
additional primary attributes that may not be persisted effectively.

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

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

=head2 oligo_type

 Title    : oligo_type
 Usage    : $obj->oligo_type
 Function : get/set the oligo type
 Returns  : the oligo type (forward_primer, reverse_primer, internal_oligo)
 Args     : optional string
 Note     : simple alias for primary_tag

=cut

=head2 rank

 Title    : rank
 Usage    : $obj->rank
 Function : get/set the rank
 Returns  : rank
 Args     : optional string

=cut

=head2 validate_seq

 Title    : validate_seq
 Usage    : $obj->validate_seq
 Function : Checks the calculated primer sequence against the actual sequence
            being analysed.
 Returns  : True (1) if validated, False (0) and a warning otherwise
 Args     : none

=cut

=head2 melting_temp

 Title    : melting_temp
 Usage    : $obj->melting_temp
 Function : returns the Tm calculated for the primer via Primer3
 Returns  : float
 Args     : optional Tm (possibly calculated via other means)

=cut

=head2 gc_content

 Title    : gc
 Usage    : $obj->gc
 Function : returns the GC content calculated for the primer via Primer3
 Returns  : float (percent)
 Args     : optional GC content (possibly calculated via other means)

=cut

=head2 run_description

 Title    : run_description
 Usage    : $obj->run_description
 Function : returns the run description for this primer (via Primer3)
 Returns  : string
 Args     : optional description

=cut

=head2 AUTOLOAD

 Title   : AUTOLOAD
 Function: Used to access tags (properties) of the primer as given by primer3
           Check the primer3 documentation for available output tags (in lower case)
           For example, the priemr3 docs describe a tag
           PRIMER_LEFT_X_SELF_ANY
           that means you can do
           my $selfbinding = $primer->self_any;
 Usage   : my $hp = $primer->hairpin_th;
           my $stab = $primer->end_stability;
           etc.
           see primer3 doc for available output tags
 Args    : none
 Returns : value for the requested tag, if tag is present, otherwise undef

=cut
