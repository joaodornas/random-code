#!/usr/bin/perl

$count = 0;

@dirs =<*_all>;
foreach $target (@dirs) {

    # print "\n------\n";
    # print "Targetdir: $target\n";

    @ROI=<$target/ROIs_*>;
    
    #    if (@ROI) {
        
        #    } else {
        
        #        print @ROI;
        #        # count++;
    
        #    }
    
    print @ROI;
    
}

# print $count
