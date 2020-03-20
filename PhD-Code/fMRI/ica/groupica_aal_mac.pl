#!/usr/bin/perl
                                                                              
use File::Basename;

$dataname="aal2std";

$aaldir="/Volumes/Dropbox/_TOOLBOX/Tim-Morten/aal2std_masks";


#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
  
    print ("\n");
  
    system ("melodic -i high-low-rest-run-1.txt -o high-low-rest-run-1-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}

foreach $i (1 .. 90){
    
    print ("\n");
    
    system ("melodic -i high-low-rest-run-2.txt -o high-low-rest-run-2-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}

exit(0);

