use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Tools::Primer3Redux',
    license             => 'perl',
    dist_author         => q{Chris Fields <cjfields@bioperl.org>},
    ##{sprintf("dist_version     => %s,",$version)##}
    build_requires => {
        'Bio::Root::Version'   => '1.006.900',
        'Test::Most'           => 0,
    },
    recursive_test_files     => 1,
    add_to_cleanup      => [ 'Bio-Tools-Primer3Redux-*' ],
);

$builder->create_build_script();
