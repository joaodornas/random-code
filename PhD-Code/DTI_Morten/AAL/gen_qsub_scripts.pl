#!/usr/bin/perl
# script to generate qsub shell scripts

use Cwd qw();
my $path = Cwd::cwd();
print "$path\n";
 
$dir=$path;
$script='fdt_script.sh';

for ($t=1; $t<91; $t++) {
    # create dir
    $regiondir=$t."_all";
    system ("mkdir $regiondir");

    #create script
    $shfile=$dir."/".$regiondir."/".$script;
    print $shfile."\n";
    open ($sh, '>', $shfile) or die "Could not create file '$shfile' $!";
    print $sh "#!/bin/sh \n";
    print $sh "cd $dir \n";
    print $sh "/usr/local/fsl/bin/probtrackx --mode=seedmask -x $dir/aal_im/aal2nodif_$t.nii.gz  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s $dir/merged -m $dir/nodif_brain_mask --dir=$dir/$regiondir --targetmasks=$dir/target_list.txt --os2t \n";
    close $sh;
    system ("chmod +x $shfile");

};
