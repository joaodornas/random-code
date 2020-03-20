#!/usr/bin/perl

# This scripts calculates the size of each ROI 

use File::Basename;

$aaldir="aal_im";

$matrixfile=">ROIs_vector.txt";
open PFILE, $matrixfile or die $!;

print "calculating size of each seed ROI\n";

foreach $i (1 .. 90){
  $target=$aaldir."/aal2nodif_".$i;
  $vol=`fslstats $target -V |awk '{print \$1}'`;

    printf PFILE $vol;
}
close (PFILE);
exit(0);



