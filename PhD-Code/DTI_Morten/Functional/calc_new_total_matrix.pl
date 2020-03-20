#!/usr/bin/perl

# NOTE: This script has to be preceded by running generate_matrix.pl
# 
# Script calculates one vector with sizes of all ROIs
# and one big connectivity matrix for ROI t_n
#  Matrix entry s_n_t_n is meanconn*ROI_no_voxels
#  where 
#    meanconn is the mean value of streamlines in voxels in ROI
#    t_n that project to s_n (max streamline in each voxel is 5000)
#    ROI_no_voxels is the number of voxels in ROI that projects from
#    s_n to t_n
#
#  Note that matrix is not symmetrical

# First calculate vector with sizes of seed ROIs
# system("calculate_vector.pl");

# Then calculate conn matrix by calculating each line and concatenating them
#  ie first use loop to calculate ROI_xx_connmatrix.txt for all ROIs
$total=0;
@dirs =<*_all>;
foreach $target (@dirs) {

    print "\n------\n";
    print "Targetdir: $target\n";

    @ROI=<$target/ROIs_*>;

    print @ROI;

    if (@ROI != ()) {
      print "OK : "; 
#      print @ROI;
      print "\n";
    }else {
      print "Not found: Calculating... \n"; 
      # calculate matrix
        #system ("cd $target; calculate_new_matrix.pl");
        system ("./calculate_new_matrix.pl");
    };

};

# Then concatenating each of the lines to one big matrix
@catfiles=();

$newfile=">newline";
open NFILE, $newfile or die $!;
printf NFILE "\n";
close (NFILE);

$matrixfile=">TOTAL_connmatrix.txt";
open PFILE, $matrixfile or die $!;

foreach $t (1..758) {

    $target=$t."_all";

    @ROI=<$target/ROIs_*>;

    push(@ROI, "newline");
    push(@catfiles, @ROI);

};

print "cat @catfiles $matrixfile";
print "\n";
system( "cat @catfiles $matrixfile");
