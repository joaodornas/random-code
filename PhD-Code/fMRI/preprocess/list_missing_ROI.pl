#!/usr/bin/perl

@dirs =<*_all>;
foreach $target (@dirs) {

    # print "\n------\n";
    # print "Targetdir: $target\n";

    #  @ROI=<$target/ROIs_*>;
    @ROI=<$target/fdt_paths*>;

    # print @ROI;
    
    print @ROI;
        
    # print @ROI;
};
