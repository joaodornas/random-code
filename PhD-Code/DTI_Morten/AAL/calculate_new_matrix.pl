#!/usr/bin/perl

$total=0;

@files =<seeds*>;
foreach $target (@files) {

    print "\n------\n";
    print "Targetseed: $target\n";

    @ii= ($target =~ m/(\d+)/g);
#    print join (",", @ii);
    $i=$ii[1];
#    print "Targetseed no: $i\n";
    $TargNum=$i;

    $vol=`fslstats $target -V |awk '{print \$1}'`;

    print "No of non-zero values: $vol";

    $mean=`fslstats $target -M`;
    print "Mean of non-zero values: $mean";

    $NumSamples=$mean*$vol;

    if ($mean>4999) {

	# this is seed roi with largest total no of connections
	$roi=$TargNum;
	print "This is the Seed ROI: $TargNum, ";
	print "Total no of samples: $NumSamples \n";
	$total=$NumSamples;

#print data $TargNum=$sum>>data

    } else {

	# this is a target ROI
	print "Target ROI: $TargNum, ";
	print "No of samples: $NumSamples";

    }
    @nsamples[$i-1]=$NumSamples;
}
print "\n";

$matrixfile=">ROIs_".$roi."_connmatrix.txt";
open PFILE, $matrixfile or die $!;

$t=1;
for $val (@nsamples) {

    # output format: 90 values separated by tabs
    $prob=$val;
    printf PFILE "%6.5f", $prob;
    print PFILE "\t";
    $t++;
};
close (PFILE);


