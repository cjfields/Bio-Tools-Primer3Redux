package Bio::Tools::Primer3Redux::PrimerPair;

use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SeqFeature::Generic);

sub left_primer {
    shift->forward_primer(@_);
}

sub forward_primer {
    my ($self, $primer) = @_;
    if ($primer) {
        $self->throw("Not a Primer object") unless $primer->isa('Bio:::Tools::Primer3Redux::Primer');
        $self->add_SeqFeature($primer, 'EXPAND');
    }
    my ($for) = grep {$_->primary_tag eq 'forward_primer'} $self->get_SeqFeatures;
    return $for;
}

sub right_primer { shift->reverse_primer(@_)}

sub reverse_primer {
    my ($self, $primer) = @_;
    if ($primer) {
        $self->throw("Not a Primer object") unless $primer->isa('Bio:::Tools::Primer3Redux::Primer');
        $self->add_SeqFeature($primer, 'EXPAND');
    }
    my ($rev) = grep {$_->primary_tag eq 'reverse_primer'} $self->get_SeqFeatures;
    return $rev;
}

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

__END__

# ABSTRACT: Simple Decorator of a Bio::SeqFeature::Generic with convenience
# methods for retrieving left and right primers, internal oligos, and any
# amplicon-related information
# AUTHOR:   Chris Fields <cjfields@cpan.org>
# OWNER:    2006-2016 Chris Fields
# LICENSE:  Perl_5
# NAME: Bio::Tools::Primer3Redux::PrimerPair - Simple Decorator of a Bio::SeqFeature::Generic with convenience
# DESCRIPTION:

Bio::Tools::Primer3Redux::PrimerPair acts as a simple SeqFeature that bundles
primer pair data together into one object.  This object can be used to retrieve
the amplicon sequence, the forward/reversion (left/right) primers, and any
internal oligos.  Furthermore, any primer information relative to the product
is included as SeqFeature tags.

# CONTRIBUTORS:

Additional contributors names and emails here

=method left_primer

 Title    : left_primer
 Note     : Alias of forward_primer()

=method forward_primer

 Title    : forward_primer
 Usage    : $obj->forward_primer
 Function : returns the forward (left) primer
 Returns  : Bio::Tools::Primer3Redux::Primer
 Args     : Optional Bio::Tools::Primer3Redux::Primer

=method right_primer

 Title    : right_primer
 Note     : alias of reverse_primer()

=method reverse_primer

 Title    : reverse_primer
 Usage    : $obj->reverse_primer
 Function : returns the reverse (right) primer
 Returns  : Bio::Tools::Primer3Redux::Primer
 Args     : Optional Bio::Tools::Primer3Redux::Primer

=method internal_oligo

 Title    : internal_oligo
 Usage    : $obj->internal_oligo
 Function : returns the internal oligo (if present)
 Returns  : Bio::Tools::Primer3Redux::Primer
 Args     : Optional Bio::Tools::Primer3Redux::Primer

