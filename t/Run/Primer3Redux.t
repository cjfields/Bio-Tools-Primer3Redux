#-*-Perl-*-
## $Id: Primer3.t 15574 2009-02-25 13:49:22Z cjfields $

# tests for Bio::Tools::Run::Primer3Rerdux
# originally written for Bio::Tools::Run::Primer3 by Rob Edwards

use strict;
use warnings;
use Data::Dumper;

BEGIN {
    use Bio::Root::Test;
    use_ok('Bio::Tools::Run::Primer3Redux');
}

use Bio::SeqIO;

my ($seqio, $seq, $primer3, $args, $results, $num_results);
$seqio = Bio::SeqIO->new(-file => test_input_file('Primer3.fa'));
$seq = $seqio->next_seq;

ok $primer3 = Bio::Tools::Run::Primer3Redux->new();

SKIP: {
    test_skip(-requires_executable => $primer3,
              -tests => 6);
    like( $primer3->program_name, qr/primer3/, 'program_name');

    # note: these are v2 parameters
    $primer3->set_parameters(
        'PRIMER_TASK'               => 'pick_pcr_primers',
        'PRIMER_SALT_CORRECTIONS'   => 1,              # optimizations
        'PRIMER_TM_FORMULA'         => 1,
    
        'PRIMER_PRODUCT_SIZE_RANGE' => '100-250',      # size
        'PRIMER_EXPLAIN_FLAG'       => 1,              # if there are errors...
    );
    
    # initial run
    ok my $parser = $primer3->run($seq);
    while (my $result = $parser->next_result) {
        isa_ok($result, 'Bio::Tools::Primer3Redux::Result');
        is($result->num_primer_pairs,5);
        my $pair = $result->next_primer_pair;
        isa_ok($pair, 'Bio::Tools::Primer3Redux::PrimerPair');
        isa_ok($pair, 'Bio::SeqFeature::Generic');

        my ($fp, $rp) = ($pair->forward_primer, $pair->reverse_primer);
        
        # can't really do exact checks here, but we can certainly check
        # various things about these...
        isa_ok($fp, 'Bio::Tools::Primer3Redux::Primer');
        isa_ok($fp, 'Bio::SeqFeature::Generic');
        cmp_ok(length($fp->seq->seq), '>', 18);
        cmp_ok(length($rp->seq->seq), '>', 18);
    }
    
    $primer3->set_parameters(
        'PRIMER_MAX_POLY_X'         => 3,   # no runs of more than 2 of same nuc
        'PRIMER_MIN_TM'             => 55,  # Tm window
        'PRIMER_MAX_TM'             => 65,
        
        'PRIMER_PRODUCT_MIN_TM'     => 75,  # product Tm
        'PRIMER_PRODUCT_OPT_TM'     => 90,
        'PRIMER_PRODUCT_MAX_TM'     => 95,
    );
    
    # Too strict!
    $parser = $primer3->run($seq);
    while (my $result = $parser->next_result) {
        isa_ok($result, 'Bio::Tools::Primer3Redux::Result');
        is($result->num_primer_pairs,0);
        
        my $pair = $result->next_primer_pair;
        ok(!defined($pair));
    }
}

done_testing();

unlink ('mlc') if -e 'mlc';
