#!/usr/bin/perl
                                                                              
use File::Basename;

$dataname="aal2std";

$aaldir="/Volumes/Dropbox/_TOOLBOX/Tim-Morten/aal2std_masks";


#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
  
    print ("\n");
  
    system ("melodic -i ../custom/filtered_func_data_mcf_unwarp2standard-clean-voxel-res -o $dataname\_$i -m $aaldir/$dataname\_$i --nobet --no_mm --report --tr=2 --Oall -v");
    
}

exit(0);

