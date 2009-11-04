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
   dist_version => '0.01',
   meta_merge => {
                  resources => {
                                repository => '/jac_sw/gitroot/perl-JCMT-MapArea.git',
                               },
                 },
   requires => {
                'perl' => '5.8.0',
                'Astro::Coords' => 0,
                'Astro::SLA' => 0,
                'Starlink::AST' => 0,
               },
   build_requires => {
                      'Test::More' => 0,
                     },
  );

$build->create_build_script;
