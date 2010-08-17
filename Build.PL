use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Tools::Primer3Redux',
    license             => 'perl',
    dist_author         => q{Chris Fields <cjfields@bioperl.org>},
    dist_version_from   => 'lib/Bio/Tools/Primer3Redux.pm',
    build_requires => {
        'Bio::Root::Version'   => '1.006.000',
        'Test::More'        => 0,
    },
    recursive_test_files     => 1,
    add_to_cleanup      => [ 'Bio-Tools-Primer3Redux-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
