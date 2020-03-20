#!/usr/bin/perl

use File::Basename;

$dataname="aal2std";

$aaldir="/home/Dropbox/_TOOLBOX/Tim-Morten/aal2std_masks";


#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
  
    print ("\n");
  
    system ("/usr/local/fsl/bin/melodic -i track-passive-trials-rest-run-1.txt -o tk-pe-tl-rt-r1-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
    
    print ("\n");
    
    system ("/usr/local/fsl/bin/melodic -i track-passive-trials-rest-run-2.txt -o tk-pe-tl-rt-r2-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
    
    print ("\n");
    
    system ("/usr/local/fsl/bin/melodic -i track-passive-trials-rest-run-3.txt -o tk-pe-tl-rt-r3-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
    
    print ("\n");
    
    system ("/usr/local/fsl/bin/melodic -i track-passive-trials-rest-run-4.txt -o tk-pe-tl-rt-r4-$dataname\_$i -m $aaldir/$dataname\_$i -a concat --nobet --no_mm --report --tr=2 --Oall -v");
    
}


exit(0);

