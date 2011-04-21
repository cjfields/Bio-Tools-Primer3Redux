# BioPerl module for Bio::Tools::Primer3Redux
#
# Copyright Chad Matsalla, Rob Edwards, Chris Fields
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer3Redux - BioPerl-based tools for parsing Primer3 output

=head1 SYNOPSIS

 # parse primer3 output to get some data
 # this is also called from Bio::Tools::Run::Primer3Redux
 use Bio::Tools::Primer3Redux;

 # read a primer3 output file
 my $p3 = Bio::Tools::Primer3Redux->new(-file=>"data/primer3_output.txt");
 
 # iterate through each result in the file
 while (my $result = $p3->next_result) {
    # iterate through each primer pair
    while (my $pair = $p3->next_primer_pair) {
        # do stuff with primer pairs...
    }
 }

=head1 DESCRIPTION

Bio::Tools::Primer3Redux creates the input files needed to design primers using
primer3 and provides mechanisms to access data in the primer3 output files.

This module provides a bioperl interface to the program primer3. See
http://www-genome.wi.mit.edu/genome_software/other/primer3.html
for details and to download the software.

This module is a reimplementation of the original Bio::Tools::Primer3 module by
Rob Edwards, itself largely based on one written by Chad Matsalla
(bioinformatics1@dieselwurks.com). Due to significant API changes between this
and the previous Primer3 module, this has been redubbed
Bio::Tools::Primer3Redux.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR -

  Chris Fields

=head1 CONTRIBUTORS

  Rob Edwards

  redwards@utmem.edu

  Brian Osborne bosborne at alum.mit.edu

  Based heavily on work of Chad Matsalla

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Primer3Redux;

our $VERSION = '0.01';

use strict;
use warnings;

use Bio::Tools::Primer3Redux::Result;

use base qw(Bio::Root::IO Bio::AnalysisParserI);

=head2 new

  Title   : new()
  Usage   : my $primer3 = Bio::Tools::Primer3->new(-file=>$file)
            to read a primer3 output file.
  Function: Parse primer3 output
  Returns : Does not return anything. If called with a filename will
            allow you to retrieve the results
  Args    : -file (optional) file of primer3 results to parse -verbose
            (optional) set verbose output
  Notes   :
  Status  : stable

=cut

=head2 next_result

  Title    : next_result
  Usage    : $obj->next_result
  Function : Returns the next Bio::Tools::Primer3Redux::Result instance (if one)
  Returns  : Bio::Tools::Primer3Redux::Result || undef
  Args     : none 
  Status   : stable
 
=cut

{

sub next_result {
	my $self = shift;
    
    $self->start_document;
    
    while (my $line = $self->_readline) {
        last if index($line, '=') == 0;
        chomp $line;
        my ($tag, $data) = split('=', $line, 2 );
        if ($tag =~ /^PRIMER_(LEFT|RIGHT|INTERNAL_OLIGO|INTERNAL|PAIR|PRODUCT)(?:(?:_(\d+))?_(.*))?/xmso) {
            my ($type, $rank, $primer_tag) = ($1, $2, $3);
            if (!defined $rank && defined $primer_tag && $primer_tag =~ /(?:(\w+)_)?(\d+)$/) {
                ($primer_tag, $rank) = ($1, $2);
            }
            $rank ||= 0;
            $type = 'INTERNAL' if $type eq 'INTERNAL_OLIGO';
            # indicates location information
            $primer_tag ||= 'LOCATION';
            if ($primer_tag eq 'EXPLAIN' || $primer_tag eq 'NUM_RETURNED') {
                $self->{persistent}->{$type}->{lc $primer_tag} = $data;
                next;
            }
            # v1 -> v2 change
            if ($type eq 'PRODUCT') {
                $type = 'PAIR';
                $primer_tag = 'PRODUCT_SIZE';
            }
            $self->{features}->{$rank}->{$type}->{lc $primer_tag} = $data;
        } elsif ($tag =~ /^(?:PRIMER_)?SEQUENCE(?:_(?:ID|TEMPLATE))?$/ )  {
            $self->{sequence}->{$tag} = $data;
        } else{ # anything else 
            $self->{run_parameters}->{$tag} = $data;
        } 
    }
    
    my $doc = $self->end_document;
    
    return $doc;
}

}

=head2 start_document

 Title    : start_document
 Usage    : $obj->start_document
 Function : initialize Primer3 output parsing
 Returns  : none
 Args     : none

=cut

sub start_document {
    my $self = shift;
    for my $data (qw(sequence features persistent run_parameters)) {
        $self->{$data} = undef;
    }
}

=head2 end_document

 Title    : end_document
 Usage    : $obj->end_document
 Function : clean up after parsing of Primer3 document is complete
 Returns  : L<Bio::Tools::Primer3Redux::Result>
 Args     : none

=cut

sub end_document {
    my $self = shift;
    my $result;
    if (defined $self->{sequence} || defined $self->{features}) {
        $result = Bio::Tools::Primer3Redux::Result->new();
        
        # data is created on the fly within the result
        $result->_initialize(-seq             => $self->{sequence},
                             -features        => $self->{features},
                             -persistent      => $self->{persistent},
                             -parameters      => $self->{run_parameters});
    }
    return $result;
}

1;

__END__

