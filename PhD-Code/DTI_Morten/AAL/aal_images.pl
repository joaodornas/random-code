#!/usr/bin/perl
$usage="# usage: aal_images.pl aal_in_dtispace
#  generates 116 nifti masks from aal (in dtispace) in aal_im dir
#  and generates a target_list.txt
#  But 
#   * Make Flirt transform from nodif_brain to std via T1
#   * Use InvertXFM to invert this nodif2std transform
#   * Use ApplyXFM using the std2nodif transform for aal into DTI space
#     using nearest neighbour interpolation.\n";
die $usage if @ARGV < 1;

                                                                              
use File::Basename;

($dataname, $path, $suffix)=fileparse(@ARGV[0], ".nii.gz");

$aaldir="aal_im";
system("mkdir $aaldir");
open (targetlist, ">>target_list.txt");

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 90){
  print ("fslmaths $dataname -thr $i -uthr $i $aaldir/$dataname\_$i");
  print ("\n");
  print (targetlist "$aaldir/$dataname\_$i.nii.gz\n");
  system ("fslmaths $dataname -thr $i -uthr $i $aaldir/$dataname\_$i");
}
close (targetlist);
exit(0);

