use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Tools::Primer3Redux',
    license             => 'perl',
    dist_author         => q{Chris Fields <cjfields@bioperl.org>},
    test_requires   => {
        'Test::Exception'   => 0,
        'Test::Warn'        => 0
        },
    build_requires => {
        'Bio::Root::Version'   => '1.006.000',
        'Test::More'        => 0,
    },
    recursive_test_files     => 1,
    add_to_cleanup      => [ 'Bio-Tools-Primer3Redux-*' ],
);

$builder->create_build_script();
