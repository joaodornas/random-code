#!/usr/bin/perl

use File::Basename;

$dataname="aal2std";

$aaldir="/Volumes/dropbox/_TOOLBOX/Tim-Morten/aal2std_masks";

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
    
    print ("\n");
    
    system ("melodic -i track-passive-trials-rest-run-4.txt -o tk-pe-tl-rt-r4-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}


exit(0);

