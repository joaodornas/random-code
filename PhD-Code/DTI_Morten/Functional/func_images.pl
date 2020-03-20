#!/usr/bin/perl
$usage="# usage: func_images.pl func_in_dtispace
#  generates 758 nifti masks from Functinal Cluster Parcellation (in dtispace) in func_im dir
#  and generates a target_list.txt
#  But 
#   * Make Flirt transform from nodif_brain to std via T1
#   * Use InvertXFM to invert this nodif2std transform
#   * Use ApplyXFM using the std2nodif transform for aal into DTI space
#     using nearest neighbour interpolation.\n";
die $usage if @ARGV < 1;

                                                                              
use File::Basename;

($dataname, $path, $suffix)=fileparse(@ARGV[0], ".nii.gz");

$funcdir="func_im";
system("mkdir $funcdir");
open (targetlist, ">>target_list.txt");

#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 758){
  print ("fslmaths $dataname -thr $i -uthr $i $funcdir/$dataname\_$i");
  print ("\n");
  print (targetlist "$funcdir/$dataname\_$i.nii.gz\n");
  system ("fslmaths $dataname -thr $i -uthr $i $funcdir/$dataname\_$i");
}
close (targetlist);
exit(0);

