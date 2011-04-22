#-*-Perl-*-
## $Id: Primer3.t 15574 2009-02-25 13:49:22Z cjfields $

# tests for Bio::Tools::Run::Primer3Redux
# originally written for Bio::Tools::Run::Primer3 by Rob Edwards

use strict;
use warnings;
use Data::Dumper;
use Bio::Root::Test;
use Bio::Tools::Run::Primer3Redux;
use Bio::SeqIO;

my ($seqio, $seq, $primer3, $args, $results, $num_results);
$seqio = Bio::SeqIO->new(-file => test_input_file('Primer3.fa'));
$seq = $seqio->next_seq;

ok $primer3 = Bio::Tools::Run::Primer3Redux->new();

my $verbose = $ENV{BIOPERLDEBUG} || 0;

my %params = (
        'PRIMER_TASK'               => 'pick_pcr_primers',
        'PRIMER_SALT_CORRECTIONS'   => 1,              # optimizations
        'PRIMER_TM_FORMULA'         => 1,
        'PRIMER_PRODUCT_SIZE_RANGE' => '100-250',      # size
        'PRIMER_EXPLAIN_FLAG'       => 1,              # if there are errors...
);

SKIP: {
    test_skip(-requires_executable => $primer3);
    like( $primer3->program_name, qr/primer3/, 'program_name');

    # note: these are v2 parameters
    $primer3->set_parameters(%params);
        
    # initial run
    ok my $parser = $primer3->run($seq);
    while (my $result = $parser->next_result) {
        isa_ok($result, 'Bio::Tools::Primer3Redux::Result');
        is($result->num_primer_pairs,5);
        
        # should be 5 unique pairs, so marshal primer locations as a key for
        # uniqueness
        my %pair_data;
        
        my $ct = 0;
        while (my $pair = $result->next_primer_pair) {
            $ct++;
            my ($fp, $rp) = ($pair->forward_primer, $pair->reverse_primer);
            my $plen = $pair->length;
            ok($plen >= 100 && $plen <=250, "product size correct: 100bp <= $plen <= 250bp");
            
            if ($ct == 1) {
                isa_ok($pair, 'Bio::Tools::Primer3Redux::PrimerPair');
                isa_ok($pair, 'Bio::SeqFeatureI');
                
                isa_ok($fp, 'Bio::Tools::Primer3Redux::Primer');
                isa_ok($fp, 'Bio::SeqFeatureI');
                
                isa_ok($rp, 'Bio::Tools::Primer3Redux::Primer');
                
                cmp_ok($fp->seq->length, '>', 18, "Length:".$fp->seq->length);
                cmp_ok($fp->gc_content, '>', 40);
                cmp_ok($fp->melting_temp, '>', 50);
                is($fp->oligo_type, 'forward_primer');
                is($fp->rank, 0);
                is($fp->run_description, 'considered 5, ok 5');
                
                cmp_ok($rp->seq->length, '>', 18);
                cmp_ok($rp->gc_content, '>', 40);
                cmp_ok($rp->melting_temp, '>', 50);
                is($rp->oligo_type, 'reverse_primer');
                is($rp->rank, 0);
                is($rp->run_description, 'considered 5, ok 5');
            }
            
            is($pair_data{sprintf("%d-%d:%d-%d",$fp->start,$fp->end,
                                  $rp->start, $rp->end)}++,
               0);
        }
        
        is(scalar(keys %pair_data), 5);
        
        # get the processed sequence, which should have the primers as SFs
        my $ps = $result->get_processed_seq;
        isa_ok($ps, 'Bio::Seq');
        is($ps->get_SeqFeatures, 15);
    }
    
    # Too strict! No results
    $primer3 = Bio::Tools::Run::Primer3Redux->new();
    $primer3->set_parameters(%params,
        'PRIMER_MAX_POLY_X'         => 3,   # no runs of more than 2 of same nuc
        'PRIMER_MIN_TM'             => 55,  # Tm window
        'PRIMER_MAX_TM'             => 65,
        
        'PRIMER_PRODUCT_MIN_TM'     => 75,  # product Tm
        'PRIMER_PRODUCT_OPT_TM'     => 90,
        'PRIMER_PRODUCT_MAX_TM'     => 95
        );
    
    $parser = $primer3->run($seq);
    while (my $result = $parser->next_result) {
        isa_ok($result, 'Bio::Tools::Primer3Redux::Result');
        is($result->num_primer_pairs,0);
        
        my $pair = $result->next_primer_pair;
        ok(!defined($pair));
        
        # sequence has no features
        my $ps = $result->get_processed_seq;
        isa_ok($ps, 'Bio::Seq');
        is($ps->get_SeqFeatures, 0);
    }
}

#unlink ('mlc') if -e 'mlc';

done_testing();
