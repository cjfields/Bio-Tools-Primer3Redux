# ABSTRACT: Create input for and work with the output from the program primer3
# AUTHOR:   Chris Fields <cjfields@cpan.org>
# OWNER:    2006-2016 Chris Fields
# LICENSE:  Perl_5

# BioPerl module for Bio::Tools::Run::Primer3Redux
#
# Copyright Chad Matsalla, Rob Edwards, Chris Fields
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 SYNOPSIS

  # design some primers.
  # the output will be put into temp.out
  use Bio::Tools::Run::Primer3Redux;
  use Bio::SeqIO;

  my $seqio = Bio::SeqIO->new(-file=>'data/dna1.fa');
  my $seq = $seqio->next_seq;

  my $primer3 = Bio::Tools::Run::Primer3Redux->new(-outfile => "temp.out",
                                            -path => "/usr/bin/primer3_core");

  # or after the fact you can change the program_name
  $primer3->program_name('my_superfast_primer3');

  unless ($primer3->executable) {
    print STDERR "primer3 can not be found. Is it installed?\n";
    exit(-1)
  }

  # set the maximum and minimum Tm of the primer
  $primer3->set_parameters('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>90);

  # Design the primers. This runs primer3 and returns a
  # Bio::Tools::Primer3::result object with the results
  # Primer3 can run in several modes (see explanation for
  # 'PRIMER_TASK' in the primer3 doccumentation). To run a task,
  # either call it by its PRIMER_TASK name as in these examples:
  $pcr_primer_results = $primer3->pick_pcr_primers($seq);
  $pcr_and_hyb_results = $primer3->pick_pcr_primers_and_hyb_probe( $seq );
  $check_results = $primer3->check_primers();

  # Alternatively, explicitly set the PRIMER_TASK parameter and
  # use the generic 'run' method (this is mainly here for backwards
  # compatibility) :
  $primer3->PRIMER_TASK( 'pick_left_only' );
  $result = $primer3->run( $seq );

  # If no task is set and the 'run' method is called, primer3 will default to
  # pick pcr primers.

  # see the Bio::Tools::Primer3Redux POD for
  # things that you can get from this. For example:

  print "There were ", $results->num_primer_pairs, " primer pairs\n";

=head1 DESCRIPTION

Bio::Tools::Run::Primer3Redux creates the input files needed to design primers
using primer3 and provides mechanisms to access data in the primer3
output files.

This module a refactoring of the original BioPerl primer3 tools, themselves a
refactoring of the original Primer3 module written by Rob Edwards. See
http://primer3.sourceforge.net for details and to download the software. This
module should work for primer3 release 1 and above but is not guaranteed to work
with earlier versions.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.html             - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 CONTRIBUTORS

Rob Edwards redwards@utmem.edu
Chad Matsalla bioinformatics1@dieselwurks.com
Shawn Hoon shawnh-at-stanford.edu
Jason Stajich jason-at-bioperl.org
Brian Osborne osborne1-at-optonline.net
Chris Fields cjfields-at-bioperl-dot-org
Frank Schwach fs5-at-sanger.ac.uk

=head1 SEE ALSO

L<Bio::Tools::Primer3>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Run::Primer3Redux;

use base qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

use strict;
use Bio::Tools::Primer3Redux;
use File::Spec;
use Scalar::Util qw(blessed reftype);
use version;

my $PROGRAMNAME = 'primer3_core';
my %PARAMS;
my @P1; #Parameters common for all 0.x and 1.x.x versions
my @P110; #Parameters for the Santalucia Tm calculations in v1.1.0-v1.1.2
my @P2; #Parameters common for all 2.x.x versions
my @P200; #2 Overlap parameters that only appears in v2.0.0
my @P220; #New and changed parameters in v2.2.0 or above
my @P222; #2 more tags added in v2.2.2
my @ALLOWED_TASKS;

# 2.0 is still in alpha (3/3/10), so fallback to v1 for determining parameters
my $DEFAULT_VERSION = version->declare('1.1.4');
BEGIN {
    # $ct assigns order of parameter building
@P1 = qw(
    PRIMER_SEQUENCE_ID
    SEQUENCE
    TARGET
    EXCLUDED_REGION
    INCLUDED_REGION
    PRIMER_COMMENT
    PRIMER_DNA_CONC
    PRIMER_EXPLAIN_FLAG
    PRIMER_FILE_FLAG
    PRIMER_FIRST_BASE_INDEX
    PRIMER_GC_CLAMP
    PRIMER_LOWERCASE_MASKING
    PRIMER_INTERNAL_OLIGO_DNTP_CONC                 PRIMER_INTERNAL_OLIGO_DIVALENT_CONC
    PRIMER_INTERNAL_OLIGO_DNA_CONC                  PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION
    PRIMER_INTERNAL_OLIGO_INPUT                     PRIMER_INTERNAL_OLIGO_MAX_GC
    PRIMER_INTERNAL_OLIGO_MAX_MISHYB                PRIMER_INTERNAL_OLIGO_MAX_POLY_X
    PRIMER_INTERNAL_OLIGO_MAX_SIZE                  PRIMER_INTERNAL_OLIGO_MAX_TM
    PRIMER_INTERNAL_OLIGO_MIN_GC                    PRIMER_INTERNAL_OLIGO_MIN_QUALITY
    PRIMER_INTERNAL_OLIGO_MIN_SIZE                  PRIMER_INTERNAL_OLIGO_MIN_TM
    PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY            PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT
    PRIMER_INTERNAL_OLIGO_OPT_SIZE                  PRIMER_INTERNAL_OLIGO_OPT_TM
    PRIMER_INTERNAL_OLIGO_SALT_CONC                 PRIMER_INTERNAL_OLIGO_SELF_ANY
    PRIMER_INTERNAL_OLIGO_SELF_END
    PRIMER_IO_WT_COMPL_ANY
    PRIMER_IO_WT_COMPL_END                          PRIMER_IO_WT_END_QUAL
    PRIMER_IO_WT_GC_PERCENT_GT                      PRIMER_IO_WT_GC_PERCENT_LT
    PRIMER_IO_WT_NUM_NS                             PRIMER_IO_WT_REP_SIM
    PRIMER_IO_WT_SEQ_QUAL                           PRIMER_IO_WT_SIZE_GT
    PRIMER_IO_WT_SIZE_LT                            PRIMER_IO_WT_TM_GT
    PRIMER_IO_WT_TM_LT
    PRIMER_LEFT_INPUT                               PRIMER_RIGHT_INPUT
    PRIMER_LIBERAL_BASE
    PRIMER_MAX_DIFF_TM                              PRIMER_MAX_END_STABILITY
    PRIMER_MAX_GC                                   PRIMER_MAX_MISPRIMING
    PRIMER_MAX_POLY_X                               PRIMER_MAX_SIZE
    PRIMER_MAX_TM
    PRIMER_MIN_END_QUALITY                          PRIMER_MIN_GC
    PRIMER_MIN_QUALITY                              PRIMER_MIN_SIZE
    PRIMER_MIN_TM
    PRIMER_MISPRIMING_LIBRARY
    PRIMER_NUM_NS_ACCEPTED                          PRIMER_NUM_RETURN
    PRIMER_OPT_GC_PERCENT                           PRIMER_OPT_SIZE
    PRIMER_OPT_TM
    PRIMER_PAIR_MAX_MISPRIMING                      PRIMER_PAIR_WT_COMPL_ANY
    PRIMER_PAIR_WT_COMPL_END                        PRIMER_PAIR_WT_DIFF_TM
    PRIMER_PAIR_WT_IO_PENALTY                       PRIMER_PAIR_WT_PRODUCT_SIZE_GT
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT                  PRIMER_PAIR_WT_PRODUCT_TM_GT
    PRIMER_PAIR_WT_PRODUCT_TM_LT                    PRIMER_PAIR_WT_PR_PENALTY
    PRIMER_PAIR_WT_REP_SIM
    PRIMER_PICK_ANYWAY                              PRIMER_PICK_INTERNAL_OLIGO
    PRIMER_PRODUCT_MAX_TM                           PRIMER_PRODUCT_MIN_TM
    PRIMER_PRODUCT_OPT_SIZE                         PRIMER_PRODUCT_OPT_TM
    PRIMER_PRODUCT_SIZE_RANGE
    PRIMER_QUALITY_RANGE_MAX                        PRIMER_QUALITY_RANGE_MIN
    PRIMER_SALT_CONC
    PRIMER_SELF_ANY                                 PRIMER_SELF_END
    PRIMER_SEQUENCE_QUALITY
    PRIMER_START_CODON_POSITION
    PRIMER_TASK
    PRIMER_WT_COMPL_ANY                             PRIMER_WT_COMPL_END
    PRIMER_WT_END_QUAL                              PRIMER_WT_END_STABILITY
    PRIMER_WT_GC_PERCENT_GT                         PRIMER_WT_GC_PERCENT_LT
    PRIMER_WT_NUM_NS                                PRIMER_WT_POS_PENALTY
    PRIMER_WT_REP_SIM                               PRIMER_WT_SEQ_QUAL
    PRIMER_WT_SIZE_GT                               PRIMER_WT_SIZE_LT
    PRIMER_WT_TM_GT                                 PRIMER_WT_TM_LT
    PRIMER_WT_TEMPLATE_MISPRIMING
    PRIMER_DEFAULT_PRODUCT                          PRIMER_DEFAULT_SIZE
    PRIMER_INSIDE_PENALTY
    PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB
    PRIMER_OUTSIDE_PENALTY
    PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
    PRIMER_MAX_TEMPLATE_MISPRIMING
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING             PRIMER_PAIR_WT_TEMPLATE_MISPRIMING
);

@P110 = qw(
    PRIMER_TM_SANTALUCIA
    PRIMER_SALT_CORRECTIONS
   PRIMER_LOWERCASE_MASKING
   PRIMER_DIVALENT_CONC
   PRIMER_DNTP_CONC
);

@P2 = qw(
    SEQUENCE_EXCLUDED_REGION
    SEQUENCE_INCLUDED_REGION
    SEQUENCE_QUALITY
    SEQUENCE_FORCE_LEFT_END
    SEQUENCE_INTERNAL_EXCLUDED_REGION
    SEQUENCE_START_CODON_POSITION
    SEQUENCE_FORCE_LEFT_START
    SEQUENCE_INTERNAL_OLIGO
    SEQUENCE_TARGET
    SEQUENCE_FORCE_RIGHT_END
    SEQUENCE_PRIMER
    SEQUENCE_TEMPLATE
    SEQUENCE_FORCE_RIGHT_START
    SEQUENCE_ID
    SEQUENCE_PRIMER_REVCOMP

    PRIMER_DNA_CONC
    PRIMER_LIBERAL_BASE
    PRIMER_PAIR_WT_PR_PENALTY
    PRIMER_DNTP_CONC
    PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING
    PRIMER_EXPLAIN_FLAG
    PRIMER_LOWERCASE_MASKING
    PRIMER_PICK_ANYWAY
    PRIMER_FIRST_BASE_INDEX
    PRIMER_MAX_END_GC
    PRIMER_PICK_INTERNAL_OLIGO
    PRIMER_GC_CLAMP
    PRIMER_MAX_END_STABILITY
    PRIMER_PICK_LEFT_PRIMER
    PRIMER_INSIDE_PENALTY
    PRIMER_MAX_GC
    PRIMER_PICK_RIGHT_PRIMER
    PRIMER_INTERNAL_DNA_CONC
    PRIMER_MAX_LIBRARY_MISPRIMING
    PRIMER_INTERNAL_DNTP_CONC
    PRIMER_MAX_NS_ACCEPTED
    PRIMER_PRODUCT_MAX_TM
    PRIMER_INTERNAL_MAX_GC
    PRIMER_MAX_POLY_X
    PRIMER_PRODUCT_MIN_TM
    PRIMER_INTERNAL_MAX_LIBRARY_MISHYB
    PRIMER_MAX_SELF_ANY
    PRIMER_PRODUCT_OPT_SIZE
    PRIMER_INTERNAL_MAX_NS_ACCEPTED
    PRIMER_MAX_SELF_END
    PRIMER_PRODUCT_OPT_TM
    PRIMER_INTERNAL_MAX_POLY_X
    PRIMER_MAX_SIZE
    PRIMER_PRODUCT_SIZE_RANGE
    PRIMER_INTERNAL_MAX_SELF_ANY
    PRIMER_MAX_TEMPLATE_MISPRIMING
    PRIMER_QUALITY_RANGE_MAX
    PRIMER_INTERNAL_MAX_SELF_END
    PRIMER_MAX_TM
    PRIMER_QUALITY_RANGE_MIN
    PRIMER_INTERNAL_MAX_SIZE
    PRIMER_MIN_END_QUALITY
    PRIMER_SALT_CORRECTIONS
    PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB
    PRIMER_MIN_GC
    PRIMER_SALT_DIVALENT
    PRIMER_INTERNAL_MAX_TM
    PRIMER_MIN_QUALITY
    PRIMER_SALT_MONOVALENT
    PRIMER_INTERNAL_MIN_GC
    PRIMER_MIN_SIZE
    PRIMER_SEQUENCING_ACCURACY
    PRIMER_INTERNAL_MIN_QUALITY
    PRIMER_MIN_THREE_PRIME_DISTANCE
    PRIMER_SEQUENCING_INTERVAL
    PRIMER_INTERNAL_MIN_SIZE
    PRIMER_MIN_TM
    PRIMER_SEQUENCING_LEAD
    PRIMER_INTERNAL_MIN_TM
    PRIMER_MISPRIMING_LIBRARY
    PRIMER_SEQUENCING_SPACING
    PRIMER_INTERNAL_MISHYB_LIBRARY
    PRIMER_NUM_RETURN
    PRIMER_TASK
    PRIMER_INTERNAL_OPT_GC_PERCENT
    PRIMER_OPT_GC_PERCENT
    PRIMER_TM_FORMULA
    PRIMER_INTERNAL_OPT_SIZE
    PRIMER_OPT_SIZE
    PRIMER_WT_END_QUAL
    PRIMER_INTERNAL_OPT_TM
    PRIMER_OPT_TM
    PRIMER_WT_END_STABILITY
    PRIMER_INTERNAL_SALT_DIVALENT
    PRIMER_OUTSIDE_PENALTY
    PRIMER_WT_GC_PERCENT_GT
    PRIMER_INTERNAL_SALT_MONOVALENT
    PRIMER_PAIR_MAX_COMPL_ANY
    PRIMER_WT_GC_PERCENT_LT
    PRIMER_INTERNAL_WT_END_QUAL
    PRIMER_PAIR_MAX_COMPL_END
    PRIMER_WT_LIBRARY_MISPRIMING
    PRIMER_INTERNAL_WT_GC_PERCENT_GT
    PRIMER_PAIR_MAX_DIFF_TM
    PRIMER_WT_NUM_NS
    PRIMER_INTERNAL_WT_GC_PERCENT_LT
    PRIMER_PAIR_MAX_LIBRARY_MISPRIMING
    PRIMER_WT_POS_PENALTY
    PRIMER_INTERNAL_WT_LIBRARY_MISHYB
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
    PRIMER_WT_SELF_ANY
    PRIMER_INTERNAL_WT_NUM_NS
    PRIMER_PAIR_WT_COMPL_ANY
    PRIMER_WT_SELF_END
    PRIMER_INTERNAL_WT_SELF_ANY
    PRIMER_PAIR_WT_COMPL_END
    PRIMER_WT_SEQ_QUAL
    PRIMER_INTERNAL_WT_SELF_END
    PRIMER_PAIR_WT_DIFF_TM
    PRIMER_WT_SIZE_GT
    PRIMER_INTERNAL_WT_SEQ_QUAL
    PRIMER_PAIR_WT_IO_PENALTY
    PRIMER_WT_SIZE_LT
    PRIMER_INTERNAL_WT_SIZE_GT
    PRIMER_PAIR_WT_LIBRARY_MISPRIMING
    PRIMER_WT_TEMPLATE_MISPRIMING
    PRIMER_INTERNAL_WT_SIZE_LT
    PRIMER_PAIR_WT_PRODUCT_SIZE_GT
    PRIMER_WT_TM_GT
    PRIMER_INTERNAL_WT_TEMPLATE_MISHYB
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT
    PRIMER_WT_TM_LT
    PRIMER_INTERNAL_WT_TM_GT
    PRIMER_PAIR_WT_PRODUCT_TM_GT
    PRIMER_INTERNAL_WT_TM_LT
    PRIMER_PAIR_WT_PRODUCT_TM_LT

    P3_FILE_ID
    P3_FILE_FLAG
    P3_COMMENT
);
@P200=qw(
    SEQUENCE_PRIMER_OVERLAP_POS
    PRIMER_POS_OVERLAP_TO_END_DIST
);

@P220=qw(
    SEQUENCE_OVERLAP_JUNCTION_LIST
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
    PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION
    PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION
    PRIMER_THERMODYNAMIC_ALIGNMENT
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH
    PRIMER_MAX_SELF_ANY_TH
    PRIMER_MAX_SELF_ANY_TH
    PRIMER_INTERNAL_MAX_SELF_ANY_TH
    PRIMER_PAIR_MAX_COMPL_ANY_TH
    PRIMER_WT_SELF_ANY_TH
    PRIMER_INTERNAL_WT_SELF_ANY_TH
    PRIMER_PAIR_WT_COMPL_ANY_TH
    PRIMER_MAX_SELF_END_TH
    PRIMER_INTERNAL_MAX_SELF_END_TH
    PRIMER_PAIR_MAX_COMPL_END_TH
    PRIMER_WT_SELF_END_TH
    PRIMER_INTERNAL_WT_SELF_END_TH
    PRIMER_PAIR_WT_COMPL_END_TH
    PRIMER_MAX_HAIRPIN_TH
    PRIMER_PAIR_MAX_HAIRPIN_TH
    PRIMER_INTERNAL_MAX_HAIRPIN_TH
    PRIMER_WT_HAIRPIN_TH
    PRIMER_INTERNAL_WT_HAIRPIN_TH
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH
    PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB_TH
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH
    PRIMER_WT_TEMPLATE_MISPRIMING_TH
    PRIMER_INTERNAL_WT_TEMPLATE_MISHYB_TH
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH

);

@P222 = qw (
    PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE
    PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE
);

# These are the allowed values for PRIMER_TASK
# A run method of the same name will be created dynamically
# in _create_run_methods for each of the tasks listed here
@ALLOWED_TASKS = qw(
  pick_pcr_primers
  pick_detection_primers
  check_primers
  pick_primer_list
  pick_sequencing_primers
  pick_cloning_primers
  pick_discriminative_primers
  pick_pcr_primers_and_hyb_probe
  pick_left_only
  pick_right_only
  pick_hyb_probe_only
);
}

=head2 new()

 Title   : new()
 Usage   : my $primer3 = Bio::Tools::Run::Primer3->new(-file=>$file) to read
           a primer3 output file.
           my $primer3 = Bio::Tools::Run::Primer3->new(-seq=>sequence object)
           design primers against sequence
 Function: Start primer3 working and adds a sequence. At the moment it
           will not clear out the old sequence, but I suppose it should.
 Returns : Does not return anything. If called with a filename will allow
           you to retrieve the results
 Args    : -outfile : file name send output results to
           -path    : path to primer3 executable
           -p3_settings_file :(optional) path to the settings file. Supported only by primer3 version 2 or above.
           -verbose :(optional) boolean value to set verbose output.


=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->io->_initialize_io();

    # create the methods for the various run methods
    $self->_create_run_methods;

    my ($program, $outfile, $path, $p3_settings_file, $verbose) = $self->_rearrange(
        [qw(PROGRAM OUTFILE PATH P3_SETTINGS_FILE VERBOSE)], @args);

    $program && $self->program_name($program);

    if ($outfile) {
        $self->outfile_name($outfile);
    }
    if ($path) {
        my (undef,$path,$prog) = File::Spec->splitpath($path);
        $self->program_dir($path);
        $self->program_name($prog);
    }

    if ($verbose) {
        $self->{'verbose'}=1;
    }
    # determine the correct set of parameters to use (v1 vs v2)

    my $v =  eval {$self->executable; 1;} ? $self->version_check : $DEFAULT_VERSION;

    #apply $p3_settings_file only if version>2
    if (($p3_settings_file)&&($v>=version->declare('2.0.0'))){
        $self->p3_settings_file($p3_settings_file);
    }

    my $ct = 0;

    #%PARAMS = ($v && $v >= version->declare("v2.0.0")) ? map {$_ => $ct++} @P2 :
    #            map {$_ => $ct++} @P1;

    if ($v){
        if ($v>=version->declare("2.2.2")){%PARAMS=map{$_ => $ct++} (@P2, @P220, @P222)}
        elsif ($v>=version->declare('2.2.0')){%PARAMS=map{$_ => $ct++} (@P2, @P220)}
        elsif ($v>=version->declare('2.0.0')){%PARAMS=map{$_ => $ct++} (@P2, @P200)}
        elsif ($v>=version->declare('1.1.0')) {%PARAMS=map{$_=> $ct++} (@P1, @P110);}
        else {%PARAMS=map{$_ => $ct++} @P1;}
    }

    $self->_set_from_args(\@args,
        -methods => [sort keys %PARAMS],
        -create => 1
       );

    return $self;
}

=head2 program_name

 Title   : program_name
 Usage   : $primer3->program_name()
 Function: holds the program name
 Returns:  string
 Args    : None

=cut

sub program_name {
    my $self = shift;
    # if explicitly set, use that
    return $self->{'program_name'} = shift @_ if @_;
    # then if previously set, use that
    return $self->{'program_name'} if $self->{'program_name'};
    # run a quick check to look for programm set class attribute if found
    if (!$PROGRAMNAME) {
        for (qw(primer3 primer3_core)) {
            if ($self->io->exists_exe($_)) {
                $PROGRAMNAME = $_;
                last;
            }
        }
    }
    # don't set permanently, use global
    return $PROGRAMNAME;
}

=head2 program_dir

 Title   : program_dir
 Usage   : $primer3->program_dir($dir)
 Function: returns the program directory, which may also be obtained from ENV variable.
 Returns :  string
 Args    :

=cut

sub program_dir {
    my ($self, $dir) = @_;
    if ($dir) {
       $self->{'program_dir'}=$dir;
    }

    # we need to stop here if we know what the answer is, otherwise we can
    # never set it and then call it later
    return $self->{'program_dir'} if $self->{'program_dir'};

    if ($ENV{PRIMER3}) {
        $self->{'program_dir'} = Bio::Root::IO->catfile($ENV{PRIMER3});
    } else {
        $self->{'program_dir'} = Bio::Root::IO->catfile('usr','local','bin');
    }

    return $self->{'program_dir'}
}

=head2  version_check

 Title   : version_check
 Usage   : $v = $prog->version_check();
 Function: Determine the version number of the program
 Example :
 Returns : a version.pm object in dotted-decimal form.
 Args    : none
 Note    : This was previously known as Bio::Tools::Run::Primer3Redux::version, and renamed to avoid confusion with version.pm

=cut

sub version_check {
    my ($self) = @_;
    return unless my $exe = $self->executable;
    if (!defined $self->{'_progversion'}) {
        my $string = `$exe -about 2>&1`;
        my $v;
        if ($string =~ m{primer3\s+release\s+([\d\.]+)}) {
            $self->{'_progversion'} = version->declare($1);
        }
    }
    return $self->{'_progversion'} || undef;
}

=head2 set_parameters()

 Title   : set_parameters()
 Usage   : $primer3->set_parameters(key=>value)
 Function: Sets parameters for the input file
 Returns : Returns the number of arguments added
 Args    : See the primer3 docs.
 Notes   : To set individual parameters use the associated method:
           $primer3->PRIMER_MAX_TM(40)

=cut

sub set_parameters {
    my ($self, %args)=@_;
    # hack around _rearrange issue to deal with lack of '-'
    my ($seq) = map {
            $args{$_}
            }
        grep { uc $_ eq 'SEQ' } keys %args;
    if (defined $seq) {
        my @seqs = (ref $seq  && (reftype $seq  eq 'ARRAY')) ? @$seq : ($seq);
        for my $s (@seqs) {
            $self->throw("-seq must be a single or array reference of Bio::SeqI") unless (blessed $seq &&
                $seq->isa('Bio::SeqI'));
        }
        $self->{seq_cache} = \@seqs;
    }

    my $added_args = 0;

    # add this back in
    unless ($self->{'no_param_checks'}) {
        for my $key (sort keys %args) {
            my $method = uc $key; # consistency
            $method =~ s/^-//;    # remove possible hanging bp-like parameter prefix
            if (!$self->can($method)) {
                next if $method eq 'SEQ';
                $self->warn("Parameter $key is not a valid Primer3 parameter");
                next
            }
            $self->$method($args{$key});
            $added_args++;
            if($self->{'verbose'}){print("$key has been updated to $args{$key}.\n");};
        }
    }
    return $added_args;
}

=head2 get_parameters

 Title    : get_parameters
 Usage    : $obj->get_parameters
 Function :
 Returns  :
 Args     :

=cut

sub get_parameters {
    my $self = shift;
    my %args = map {$_->[0] => $_->[1]}
        grep { defined $_->[1] }
        map { [$_, $self->$_] } sort keys %PARAMS;
    return %args;
}

=head2 reset_parameters()

 Title   : reset_parameters()
 Usage   : $primer3->reset_parameters()
 Function: Resets all parameters to be undef
 Returns : none
 Args    : none; to reset specific targets call the specific method for that
           target (i.e. $primer3->PRIMER_MAX_TM(undef))

=cut

sub reset_parameters {
    my $self = shift;
    my %args = map {$_ => undef} sort keys %PARAMS;
    $self->_set_from_args(\%args, -methods => [sort keys %PARAMS]);
}

=head2 p3_settings_file()

 Title  : p3settingsfile()
 Usage  : $primer3->p3settingsfile($file_path);
 Function   : Getter/Setter for the Primer3 settings file.
 Returns    : A string containing the path of the named settings file.
 Args   : $file_path A valid file path to the settings file.
 Note   : This argument only works in primer3 version 2 or above.

=cut

sub p3_settings_file{
    my ($self, $file_path)=@_;
        if ($self->version_check()>=version->declare('2.0.0')){ #version check
            if ($file_path){$self->{'p3_settings_file'}=$file_path;}
            if ($self->{'p3_settings_file'}){return $self->{'p3_settings_file'};}
        }
        else {$self->warn("p3_settings_file is only available on primer3 release 2 or above.")}
}

=head2 run

 Title   : run
 Usage   : $primer3->run;
 Function: Generic run method for the primer3 program. The PRIMER_TASK
           must be defined either in the parameter set or it defaults to
           'pick_pcr_primers'.
 Returns : A Bio::Tools::Primer3 object containing the results.
           See the Bio::Tools::Primer3 documentation for those functions.
 Args    : Bio::Seq object(s) to use as SEQUENCE_TEMPLATE(s)

=cut

sub run {
   my($self, @seqs) = @_;
   $self->_do_run(\@seqs);
}


# Create a method that calls _do_run with correct parameters
# for each of the allowed run methods (tasks)
# so we can run primer3 with calls like this:
# $primer3->pick_left_only( $seq );

sub _create_run_methods {
    my $self = shift;
    foreach my $task (@ALLOWED_TASKS) {
      next if defined  *{ __PACKAGE__ . '::' . $task };
        my $code = 'my($self, @seqs) = @_; $self->_do_run(\@seqs, "'.
          $task.'");';
        my $sub = eval "sub { $code }";
        $self->throw("Compilation error while creating method for $task: $@") if $@;
        no strict 'refs';
        *{ __PACKAGE__ . '::' . $task } = $sub;
    }
}



=head2 _do_run

 Title   : _do_run
 Usage   : INTERNAL
 Function: Generate input file and run the primer3 program.
 Returns : A Bio::Tools::Primer3 object containing the results.
           See the Bio::Tools::Primer3 documentation for those functions.
 Args    : sequence or ref to array of sequences, [task]
           If 'task' given, parameter 'PRIMER_TASK' is set to the given task,
           otherwise, the 'PRIMER_TASK' parameter is used to detmine what to do.
           If that too isn't set the primer3 will default to pick_pcr_primers

=cut

sub _do_run {
    my($self, $seqs, $task) = @_;
    my @seqs = ref $seqs eq 'ARRAY' ? @$seqs : ($seqs);

    my @exec_array;
    my $executable = $self->executable;
    my $arguments = ''; #variable to hole run-time arguments
    my $out = $self->outfile_name;
    unless ($executable && -e $executable) {
        $self->throw("Executable was not found. Do not know where primer3 is!") if !$executable;
        $self->throw("$executable was not found. Do not know where primer3 is!");
        exit(-1);
    }

    push (@exec_array, $executable);

    $self->PRIMER_TASK( $task ) if $task;
    my %params = $self->get_parameters;

    my $file = $self->_generate_input_file(\%params, \@seqs);

    if($self->version_check>= version->declare('v2.0.0')){
        if ($self->p3_settings_file()){push (@exec_array, "-p3_settings_file=". $self->p3_settings_file());}
        if ($self->{'verbose'}){
            push (@exec_array, "-echo_settings_file");
        }
    }

    my $str;
    foreach (@exec_array){$str.=$_." ";} #"$executable $arguments < $file";

    my $obj = Bio::Tools::Primer3Redux->new(-verbose => $self->verbose);

    my $status;
    my @args;
    # file output
    if ($out) {
        if($self->version_check>=version->declare('v2.0.0')){
            $str.= "-output=".$out. " <". $file;
            $status=system ($str);
        }
        else {
            $str .= "< $file > $out";
        my $status = system($str);
        }
        if($status || !-e $out || -z $out ) {
            my $error = ($!) ? "$! Status: $status" : "Status: $status";
            $self->throw( "Primer3 call crashed: $error \n[command $str]\n");
            return undef;
        }
        if ($obj && ref($obj)) {
            $obj->file($out);
            @args = (-file => $out);
        }
    # fh-based (no outfile)
    } else {
        $str.=' <'. $file;
        open(my $fh,$str."|") || $self->throw("Primer3 call ($str) crashed: $?\n");
        if ($obj && ref($obj)) {
            $obj->_fh($fh);
            @args = (-fh => $fh);
        } else {
            # dump to debugging
            my $io;
            while(<$fh>) {$io .= $_;}
            close($fh);
            $self->debug($io);
            return 1;
        }
    }
    $obj->_initialize_io(@args) if $obj && ref($obj);
    return $obj;
}


sub _generate_input_file {
    # note that I write this to a temp file because we need both read
    # and write access to primer3, therefore,
    # we can't use a simple pipe.
    my ($self, $args, $seqs) = @_;
    my ($tmpfh, $tmpfile) = $self->io->tempfile();

    # this is a hack to get around interface issues and conflicts when passing
    # in raw sequence via PRIMER_SEQUENCE_ID and SEQUENCE (one can potentially
    # have both).  For now, push any explicitly set parameters on last

    my ($id_tag, $seq_tag) = $self->version_check >= version->declare ('v2.0.0') ? # v2 differs from v1
                            qw(SEQUENCE_ID SEQUENCE_TEMPLATE) :  #v2
                            qw(PRIMER_SEQUENCE_ID SEQUENCE);   #v1
    my @seqdata;
    for my $seq (@$seqs) {
        $self->throw("Arguments to run() must be a single or array ref of Bio::SeqI")
            unless (blessed $seq && $seq->isa('Bio::SeqI')) ;
        push @seqdata, {$id_tag     => $seq->id,
                        $seq_tag    => $seq->seq};
    }

    if (exists $args->{$id_tag} || exists $args->{$seq_tag}) {
        push @seqdata, {$id_tag     => $args->{$id_tag},
                        $seq_tag    => $args->{$seq_tag}};
        delete $args->{$id_tag};
        delete $args->{$seq_tag};
    }

    # generate the common BoulderIO string to be used for each sequence
    my $string = '';

    for my $param (sort keys %$args) {
        my $tmp = $self->$param;
        my @data = (ref $tmp && reftype $tmp eq 'ARRAY' ) ? @$tmp : $tmp;
        for my $d (@data) {
            $string .= "$param=$d\n";
        }
    }

    $string .= "=\n";

    # Print the parameters once for each sequence provided or, if
    # no sequence given, just write a single block without the
    # SEQUENCE_TEMPLATE but make sure we are doing a check_primers taks
    # because that's the only allowed situation where no sequence is required
    if (@seqdata){
      for my $data (@seqdata) {
          my $str = join("\n", map { "$_=".$data->{$_}} sort keys %$data)."\n$string";
          if($self->{'verbose'}){$self->debug("TRYING\n$str");}
          print $tmpfh $str;
      }
    } elsif ( $args->{PRIMER_TASK} eq 'check_primers' ){
      print $tmpfh $string;
    } else {
      $self->throw("'run' was called without a sequence but PRIMER_TASK was not ".
        "'check_primers'. All tasks except 'check_primers' must be given a sequence ".
        "to design primers on");
    }

    close($tmpfh);
    return $tmpfile;
}

1;
