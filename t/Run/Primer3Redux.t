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
$seqio=Bio::SeqIO->new(-file => test_input_file('Primer3.fa'));
$seq=$seqio->next_seq;

# sequence is passed to run() now
ok $primer3 = Bio::Tools::Run::Primer3Redux->new(); 

SKIP: {
    test_skip(-requires_executable => $primer3,
              -tests => 6);
    
    $args = $primer3->arguments;
    is($$args{'PRIMER_SEQUENCE_ID'}, "(string, optional) an id. Optional. Note must be present if PRIMER_FILE_FLAG is set");
    ##ok $primer3->add_targets('PRIMER_SEQUENCE_ID'=>'test seq');
    ok my $parser = $primer3->run($seq);
    print STDERR $parser."\n";
    while (my $result = $parser->next_result) {
        
    #    is($results->number_of_results,5);
    #    is( $results->{input_options}->{PRIMER_SEQUENCE_ID},'test seq');
    }
    like( $primer3->program_name, qr/primer3/, 'program_name');
}

done_testing();
