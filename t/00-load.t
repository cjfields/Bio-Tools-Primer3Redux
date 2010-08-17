#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Bio::Tools::Primer3Redux' ) || print "Bail out!
";
}

diag( "Testing Bio::Tools::Primer3Redux $Bio::Tools::Primer3Redux::VERSION, Perl $], $^X" );
