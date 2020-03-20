#!/usr/bin/perl

$basedir=".";
$seeddir="func_im";

$targetfile="target_list.txt";

for $i (1..758) {

$seed="FuncClusterParcels2nodif_$i.nii.gz";
$seedmask=$basedir."/".$seeddir."/".$seed;

$outbasedir=$i."_all";
$outputdir=$basedir."/".$outbasedir;
$tmasks=$basedir."/".$targetfile;
$nodifmask=$basedir."/"."nodif_brain_mask";
$merged=$basedir."/"."merged";

print ("/usr/local/fsl/bin/probtrackx --mode=seedmask -x $seedmask  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s $merged -m $nodifmask --dir=$outputdir --targetmasks=$tmasks --os2t");
print "\n";
system ("/usr/local/fsl/bin/probtrackx --mode=seedmask -x $seedmask  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s $merged -m $nodifmask --dir=$outputdir --targetmasks=$tmasks --os2t");

};
