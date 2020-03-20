#!/usr/bin/perl

# script to run qsub shell scripts
 
use Cwd qw();
my $path = Cwd::cwd();
print "$path\n";
 
$dir=$path;
$script='fdt_script.sh';

for ($t=1; $t<759; $t++) {
    # create dir
    $regiondir=$t."_all";
    print ("qsub -V -cwd -e e$t -o i$t $regiondir/$script\n");
    system ("qsub -V -cwd -e e$t -o i$t $regiondir/$script");
};
