package JCMT::MapArea;

=head1 NAME

JCMT::MapArea - Calculate requested map area for JCMT instruments.

=head1 SYNOPSIS

  use JCMT::MapArea;
  my $region = JCMT::MapArea::region( header => $header );

=head1 DESCRIPTION

This class calculates an AST Region delineating the requested map area
based on FITS headers from a JCMT observation.

=cut

use strict;
use warnings;
use Carp;

our $DEBUG = 0;
our $VERSION = '0.01';

use Astro::Coords;
use Astro::PAL qw/ :pal :constants /;
use Starlink::AST;

my %TRACK2AST = ( 'J2000' => 'FK5',
                  'B1950' => 'FK4',
                  'APP' => 'GAPPT',
                  'GAL' => 'GALACTIC' );

=head1 METHODS

=over 4

=item B<region>

Return an AST Region based on FITS headers.

  my $region = region( header => $header );

The named argument 'header' must be a hash reference, with keys being
header keys and values being header values.

=cut

sub region {
  my %args = @_;
  if( ! defined( $args{'header'} ) ) {
    return undef;
  }

  my $hdr = $args{'header'};

  # Retrieve header values.
  my $basec1   = _getval( $hdr, "BASEC1" );
  my $basec2   = _getval( $hdr, "BASEC2" );
  my $map_x    = $hdr->{ "MAP_X" };
  my $map_y    = $hdr->{ "MAP_Y" };
  my $map_pa   = $hdr->{ "MAP_PA" };
  my $map_hght = $hdr->{ "MAP_HGHT" };
  my $map_wdth = $hdr->{ "MAP_WDTH" };
  my $tracksys = uc( $hdr->{ "TRACKSYS" } );

  if( ! defined( $basec1 ) ) {
    croak "Must define BASEC1 to calculate AST Region";
  }
  if( ! defined( $basec2 ) ) {
    croak "Must define BASEC2 to calculate AST Region";
  }
  if( ! defined( $map_x ) ) {
    croak "Must define MAP_X to calculate AST Region";
  }
  if( ! defined( $map_y ) ) {
    croak "Must define MAP_Y to calculate AST Region";
  }
  if( ! defined( $map_pa ) ) {
    croak "Must define MAP_PA to calculate AST Region";
  }
  if( ! defined( $map_hght ) ) {
    croak "Must define MAP_HGHT to calculate AST Region";
  }
  if( ! defined( $map_wdth ) ) {
    croak "Must define MAP_WDTH to calculate AST Region";
  }
  if( ! defined( $tracksys ) ) {
    carp "TRACKSYS not defined; defaulting to J2000";
    $tracksys = 'J2000';
  }

  my $date_obs = _getval( $hdr, "DATE-OBS" );
  if( $tracksys eq 'APP' && ! defined( $date_obs ) ) {
    croak "When TRACKSYS is APP, DATE_OBS must be defined";
  }

  my ( $basera, $basedec );

  if ($tracksys eq 'GAL') {
    my $base = new Astro::Coords( long => $basec1,
                                  lat => $basec2,
                                  type => $tracksys,
                                  units => 'degrees' );
    ($basera, $basedec) = $base->glonglat;

  } elsif ($tracksys eq 'APP') {
    my $base = new Astro::Coords( ra1 => $basec1,
                                  ra2 => $basec1,
                                  dec1 => $basec2,
                                  dec2 => $basec2,
                                  units => 'degrees',
                                  mjd1 => 0,
                                  mjd2 => 0,
                                );
    ( $basera, $basedec ) = $base->apparent();

  } else {
    my $base = new Astro::Coords( ra => $basec1,
                                  dec => $basec2,
                                  type => $tracksys,
                                  units => 'degrees' );
    if ($tracksys eq 'B1950') {
      ($basera, $basedec) = $base->radec1950();
    } else {
      ($basera, $basedec) = $base->radec();
    }
  }

# Calculate the rotation angle in radians. Note that this position
# angle is north of east whereas that in the header is east of north.
  my $rot = -1.0 * $map_pa * DD2R;

# Calculate the X and Y positions of the map in the map's rotation
# frame.
  my $mx1 =        $map_wdth / 2;
  my $my1 =        $map_hght / 2;
  my $mx2 =        $map_wdth / 2;
  my $my2 = -1.0 * $map_hght / 2;
  my $mx3 = -1.0 * $map_wdth / 2;
  my $my3 = -1.0 * $map_hght / 2;
  my $mx4 = -1.0 * $map_wdth / 2;
  my $my4 =        $map_hght / 2;

# Convert these positions into rotated offsets.
  my ( $mrx1, $mry1 ) = _rotate( $mx1, $my1, $rot );
  my ( $mrx2, $mry2 ) = _rotate( $mx2, $my2, $rot );
  my ( $mrx3, $mry3 ) = _rotate( $mx3, $my3, $rot );
  my ( $mrx4, $mry4 ) = _rotate( $mx4, $my4, $rot );

# Calculate the rotated offset position.
  my ( $rx, $ry ) = _rotate( $map_x, $map_y, $rot );

# Calculate the corners in spherical coordinates.
  my ($c1ra, $c1dec) = palDtp2s( ( $rx + $mrx1 ) * DAS2R,
                                 ( $ry + $mry1 ) * DAS2R,
                                 $basera, $basedec );
  my ($c2ra, $c2dec) = palDtp2s( ( $rx + $mrx2 ) * DAS2R,
                                 ( $ry + $mry2 ) * DAS2R,
                                 $basera, $basedec );
  my ($c3ra, $c3dec) = palDtp2s( ( $rx + $mrx3 ) * DAS2R,
                                 ( $ry + $mry3 ) * DAS2R,
                                 $basera, $basedec );
  my ($c4ra, $c4dec) = palDtp2s( ( $rx + $mrx4 ) * DAS2R,
                                 ( $ry + $mry4 ) * DAS2R,
                                 $basera, $basedec );

# Create the region.
  my $ast_system = $TRACK2AST{ uc( $tracksys ) };
  my $params = "SYSTEM=$ast_system";
  if( $tracksys eq 'APP' ) {
    $params .= ",EPOCH=$date_obs";
  }
  my $skyFrame = new Starlink::AST::SkyFrame( $params );
  my $region = Starlink::AST::Polygon->new( $skyFrame,
                                            [ $c1ra, $c2ra, $c3ra, $c4ra ],
                                            [ $c1dec, $c2dec, $c3dec, $c4dec ],
                                            undef, "" );

# Return it.
  return $region;

}

=begin __PRIVATE__

=item B<_rotate>

Rotate an X and Y position through a given angle, returning the
rotated position.

 ( $xprime, $yprime ) = _rotate( $x, $y, $rot );

Rotation angle must be in radians.

=cut

sub _rotate {
  my $x = shift;
  my $y = shift;
  my $rot = shift;

  my $xprime = $x * cos( $rot ) - $y * sin( $rot );
  my $yprime = $x * sin( $rot ) + $y * cos( $rot );

  return ( $xprime, $yprime );
}

=item B<_getval>

Retrieve value from the header, including subheaders.

=cut

sub _getval {
  my $hdr = shift;
  my $key = shift;

  if (exists $hdr->{$key}) {
    return $hdr->{$key};
  }

  if (exists $hdr->{SUBHEADERS}) {
    for my $sh (@{$hdr->{SUBHEADERS}}) {
      if (exists $sh->{$key} && defined $sh->{$key}) {
        return $sh->{$key};
      }
    }
  }
  return;
}

=end __PRIVATE__

=back

=head1 AUTHOR

Brad Cavanagh E<lt>b.cavanagh@jach.hawaii.eduE<gt>

=head1 COPYRIGHT

Copyright 2009 Science and Technology Facilities Council. All Rights
Reserved.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place,Suite 330, Boston, MA 02111-1307,
USA

=cut

1;
