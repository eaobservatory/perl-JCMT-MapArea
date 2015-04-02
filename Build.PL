#!perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new
  (
   module_name => 'JCMT::MapArea',
   license => 'gpl',
   dist_name => 'JCMT::MapArea',
   dist_author => [
                   'Brad Cavanagh <b.cavanagh@jach.hawaii.edu>',
                  ],
   dist_version => '0.02',
   meta_merge => {
                  resources => {
                                repository => 'git://github.com/eaobservatory/perl-JCMT-MapArea.git',
                               },
                 },
   requires => {
                'perl' => '5.8.0',
                'Astro::Coords' => 0,
                'Astro::PAL' => 0,
                'Starlink::AST' => 0,
                'Astro::FITS::Header' => 0,
               },
   build_requires => {
                      'Test::More' => 0,
                     },
   script_files => [ qw/ maparea2region / ],
  );

$build->create_build_script;
