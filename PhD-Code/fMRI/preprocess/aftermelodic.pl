#!/usr/bin/perl

$dir{1}="T2-Stimulus-MOT-4-Balls-Track/Run-1-1-1/FSL/custom";
$dir{2}="teste2";


#only include cortical and subcortical areas, NOT cerebellum
foreach $i (1 .. 2){
  
  system ("mkdir $dir{$i}");
    
}


