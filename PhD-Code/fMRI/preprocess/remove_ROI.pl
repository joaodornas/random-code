#!/usr/bin/perl

@dirs =<*_all>;
foreach $target (@dirs) {

    # print "\n------\n";
    # print "Targetdir: $target\n";

    @ROI=<$target/ROIs_*>;

    # print @ROI;
    
    foreach $this_ROI (@ROI) {
        
        system ("rm $this_ROI");

    }
};
