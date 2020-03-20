#!/usr/bin/perl

for $i (0..166) {

    if ($i < 10) {
        $j = "000$i";
    }
    if (($i >= 10) && ($i < 100)) {
        $j = "00$i";
    }
    if ($i >= 100) {
        $j = "0$i";
    }
    
system ("avscale --allparams prefiltered_func_data_mcf.mat/MAT_$j > MAT_$j.txt");

};
