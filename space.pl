#!/usr/bin/env perl

=head1 NAME

B<space_molding>

=head1 DESCRIPTION

A novel space aided homology modeling approach.

=head1 SYNOPSIS

./space_molding.pl file

=head1 DEPENDENCIES

-Perl, >= 5.10

-Bio::PDB::Structure, >= 0.01 (Perl module)

-File::Spec, >= 3.33 (Perl module)

-List::Util, >= 1.25 (Perl module)

-Math::Round, >= 0.06 (Perl module)

=head1 AUTHORS

-Dimitrios - Georgios Kontopoulos <<dgkontopoulos@gmail.com>>

-Dimitrios Vlachakis <<dvlachakis@bioacademy.gr>>

-Sophia Kossida <<skossida@bioacademy.gr>>

=head1 LICENSE

This program is free software: you can redistribute it 
and/or modify it under the terms of the GNU General 
Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your 
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For more information, see
<a href="http://www.gnu.org/licenses/" style="text-decoration:none">
http://www.gnu.org/licenses/<a>.

=cut

use strict;
use warnings;

use feature qw(say);

use Bio::PDB::Structure;
use File::Spec;
use List::Util qw(min max);
use Math::Round qw(round);

@ARGV == 1 or die "Usage: $0 file\n";

#File and directory handling.#
my $input_file = $ARGV[0];

$input_file = File::Spec->rel2abs($input_file);
open my $fh, '<', $input_file or die "ERROR. $input_file cannot be opened.\n";
close $fh;

my $filename;
( undef, undef, $filename ) = File::Spec->splitpath($input_file);

my $directory;
if ( $input_file =~ /$filename$/ )
{
    $directory = $`;
}

my $prot = Bio::PDB::Structure::Molecule->new;

my ( $atoms_number, %grid );

my $xmin = 99_999;
my $ymin = 99_999;
my $zmin = 99_999;
my $xmax = -99_999;
my $ymax = -99_999;
my $zmax = -99_999;

if ( $input_file =~ /[.]pdb$/ )
{

    #Read the file.#
    $prot->read($input_file);
    $atoms_number = $prot->size;

    #For every atom...#
    for my $atom ( 0 .. $atoms_number - 1 )
    {

        #... round its coordinates...#
        my $atom1  = $prot->atom($atom);
        my $x_coor = round( $atom1->x );
        my $y_coor = round( $atom1->y );
        my $z_coor = round( $atom1->z );

        #... and put them in the grid.#
        if ( !( defined $grid{"$x_coor $y_coor $z_coor"} ) )
        {
            $grid{"$x_coor $y_coor $z_coor"} = 1;
        }

        #Get the current grid dimensions.#
        $xmin = min $xmin, $x_coor;
        $xmax = max $xmax, $x_coor;
        $ymin = min $ymin, $y_coor;
        $ymax = max $ymax, $y_coor;
        $zmin = min $zmin, $z_coor;
        $zmax = max $zmax, $z_coor;
    }
}

#Create the model.#
my $output =
  modeling(
    proximity_check( $xmin, $xmax, $ymin, $ymax, $zmin, $zmax, \%grid ) );

say "\nA model was successfully created at \"$output\"!\n";

#######################################################
# S    U    B    R    O    U    T    I    N    E    S #
#######################################################

#Select air atoms located next to one and only atom of the protein.#
sub proximity_check
{
    my ( $xmin, $xmax, $ymin, $ymax, $zmin, $zmax, $grid ) = @_;

    #Iterate through the grid elements.#
    for my $current_x ( $xmin - 1 .. $xmax + 1 )
    {
        for my $current_y ( $ymin - 1 .. $ymax + 1 )
        {
            for my $current_z ( $zmin - 1 .. $zmax + 1 )
            {
                if ( !( defined $grid{"$current_x $current_y $current_z"} ) )
                {

                    #Prepare the steps.#
                    my $temp_x_p = $current_x + 1;
                    my $temp_y_p = $current_y + 1;
                    my $temp_z_p = $current_z + 1;
                    my $temp_x_m = $current_x - 1;
                    my $temp_y_m = $current_y - 1;
                    my $temp_z_m = $current_z - 1;

                    my $sum = 0;

                    my @neighbors = (
                        "$temp_x_p $current_y $current_z",
                        "$temp_x_p $temp_y_p $current_z",
                        "$temp_x_p $temp_y_m $current_z",
                        "$temp_x_m $current_y $current_z",
                        "$temp_x_m $temp_y_p $current_z",
                        "$temp_x_m $temp_y_m $current_z",
                        "$current_x $temp_y_p $current_z",
                        "$current_x $temp_y_p $temp_z_p",
                        "$current_x $temp_y_p $temp_z_m",
                        "$current_x $temp_y_m $temp_z_m",
                        "$current_x $temp_y_m $temp_z_p",
                        "$current_x $temp_y_m $current_z",
                        "$current_x $current_y $temp_z_p",
                        "$temp_x_p $current_y $temp_z_p",
                        "$temp_x_m $current_y $temp_z_p",
                        "$temp_x_m $current_y $temp_z_m",
                        "$temp_x_p $current_y $temp_z_m",
                        "$current_x $current_y $temp_z_m",
                    );

                    for ( 0 .. $#neighbors )
                    {

                        #Check for existing protein atoms nearby.#
                        if ( $grid->{ $neighbors[$_] } )
                        {
                            $sum += $grid->{ $neighbors[$_] };
                            last if ( $sum > 1 );
                        }
                    }

                    #Select only those having one protein atom next to them.#
                    if ( $sum == 1 )
                    {
                        $grid->{"$current_x $current_y $current_z"} = 0;
                    }
                }
            }
        }
    }
    return $grid;
}

#Generate a PDB file with dummy atoms.#
sub modeling
{
    my ($grid) = @_;

    my ( $x_coor, $y_coor, $z_coor );
    my $atom_number = 1;

    if ( $filename =~ /[.]pdb/ )
    {
        $filename = $`;
    }

    #Set the output file's name and change it if it already exists.#
    my $resulting_file = $directory . $filename . '_model.pdb';

    while ( -e $resulting_file )
    {
        if ( $resulting_file =~ /([_]?\w+)[.]pdb/ )
        {
            $resulting_file = $` . $1 . '_new.pdb';
        }
    }

    open my $fh, '>', $resulting_file or die "Cannot create $resulting_file.\n";
    foreach my $key ( keys %{$grid} )
    {
        if ( $grid{$key} == 0 )
        {
            if ( $key =~ /([-]?\d+)\s+([-]?\d+)\s+([-]?\d+)/ )
            {

                #Make sure the PDB file is well formatted.#
                $x_coor = $1;
                $x_coor = sprintf '%7.3f', $x_coor;
                $y_coor = $2;
                $y_coor = sprintf '%7.3f', $y_coor;
                $z_coor = $3;
                $z_coor = sprintf '%7.3f', $z_coor;
            }
            my $atom_number_new = sprintf '%4d', $atom_number;

            #Print to the resulting file.#
            say {$fh}
"HETATM $atom_number_new Du       U         $x_coor $y_coor $z_coor  0.00";
            $atom_number++;
        }
    }
    say {$fh} 'END';
    close $fh;
    return $resulting_file;
}
