function seeds = getFunctionalSeeds_v3(network_label)

%%% coordinates in mm

if strcmp(network_label,'DAN') %%% 100
   
    seeds.network_label = 'DAN';
    
    seeds.ROI(1).label = 'r-FEF-1';
    seeds.ROI(1).x = 30.3;
    seeds.ROI(1).y = -12.8;
    seeds.ROI(1).z = 52.6;
    seeds.ROI(1).idx = 101;

    seeds.ROI(2).label = 'r-FEF-2';
    seeds.ROI(2).x = 29.0;
    seeds.ROI(2).y = -8.0;
    seeds.ROI(2).z = 56.0;
    seeds.ROI(2).idx = 102;
    
    seeds.ROI(3).label = 'r-FEF-3';
    seeds.ROI(3).x = 27.0;
    seeds.ROI(3).y = -14.0;
    seeds.ROI(3).z = 58.0;
    seeds.ROI(3).idx = 103;
    
    seeds.ROI(4).label = 'l-FEF-1';
    seeds.ROI(4).x = -26.3;
    seeds.ROI(4).y = -11.8;
    seeds.ROI(4).z = 52.7;
    seeds.ROI(4).idx = 104;

    seeds.ROI(5).label = 'l-FEF-2';
    seeds.ROI(5).x = -34.0;
    seeds.ROI(5).y = -13.0;
    seeds.ROI(5).z = 58.0;
    seeds.ROI(5).idx = 105;
    
    seeds.ROI(6).label = 'l-FEF-3';
    seeds.ROI(6).x = -25.0;
    seeds.ROI(6).y = -16.0;
    seeds.ROI(6).z = 61.0;
    seeds.ROI(6).idx = 106;

    seeds.ROI(7).label = 'r-AIPS';
    seeds.ROI(7).x = 34.0;
    seeds.ROI(7).y = -51.0;
    seeds.ROI(7).z = 49.0;
    seeds.ROI(7).idx = 107;

    seeds.ROI(8).label = 'l-AIPS';
    seeds.ROI(8).x = -42.0;
    seeds.ROI(8).y = -44.0;
    seeds.ROI(8).z = 46.0;
    seeds.ROI(8).idx = 108;

    seeds.ROI(9).label = 'r-PIPS-1';
    seeds.ROI(9).x = 23.3;
    seeds.ROI(9).y = -69.4;
    seeds.ROI(9).z = 48.6;
    seeds.ROI(9).idx = 109;
    
    seeds.ROI(10).label = 'r-PIPS-2';
    seeds.ROI(10).x = 31.0;
    seeds.ROI(10).y = -65.0;
    seeds.ROI(10).z = 57.0;
    seeds.ROI(10).idx = 110;
    
    seeds.ROI(11).label = 'r-PIPS-3';
    seeds.ROI(11).x = 19.0;
    seeds.ROI(11).y = -72.0;
    seeds.ROI(11).z = 57.0;
    seeds.ROI(11).idx = 111;
    
    seeds.ROI(12).label = 'l-PIPS-1';
    seeds.ROI(12).x = -25.3;
    seeds.ROI(12).y = -67.3;
    seeds.ROI(12).z = 47.6;
    seeds.ROI(12).idx = 112;
    
    seeds.ROI(13).label = 'l-PIPS-2';
    seeds.ROI(13).x = -29.0;
    seeds.ROI(13).y = -64.0;
    seeds.ROI(13).z = 55.0;
    seeds.ROI(13).idx = 113;
    
    seeds.ROI(14).label = 'l-PIPS-3';
    seeds.ROI(14).x = -22.0;
    seeds.ROI(14).y = -72.0;
    seeds.ROI(14).z = 51.0;
    seeds.ROI(14).idx = 114;
    
    seeds.ROI(15).label = 'r-VIPS';
    seeds.ROI(15).x = 30.3;
    seeds.ROI(15).y = -83.2;
    seeds.ROI(15).z = 13.0;
    seeds.ROI(15).idx = 115;
    
    seeds.ROI(16).label = 'l-VIPS-1';
    seeds.ROI(16).x = -24.2;
    seeds.ROI(16).y = -72.6;
    seeds.ROI(16).z = 28.8;
    seeds.ROI(16).idx = 116;
    
    seeds.ROI(17).label = 'l-VIPS-2';
    seeds.ROI(17).x = -25.0;
    seeds.ROI(17).y = -76.0;
    seeds.ROI(17).z = 26.0;
    seeds.ROI(17).idx = 117;
    
    seeds.ROI(18).label = 'r-SPL';
    seeds.ROI(18).x = 13.0;
    seeds.ROI(18).y = -72.0;
    seeds.ROI(18).z = 60.0;
    seeds.ROI(18).idx = 118;
    
    seeds.ROI(19).label = 'l-SPL';
    seeds.ROI(19).x = -14.0;
    seeds.ROI(19).y = -74.0;
    seeds.ROI(19).z = 53.0;
    seeds.ROI(19).idx = 119;
    
    seeds.ROI(20).label = 'r-MT-1';
    seeds.ROI(20).x = 42.4;
    seeds.ROI(20).y = -69.7;
    seeds.ROI(20).z = -11.2;
    seeds.ROI(20).idx = 120;
    
    seeds.ROI(21).label = 'r-MT-2';
    seeds.ROI(21).x = 51.0;
    seeds.ROI(21).y = -63.0;
    seeds.ROI(21).z = -14.0;
    seeds.ROI(21).idx = 121;
    
    seeds.ROI(22).label = 'l-MT-1';
    seeds.ROI(22).x = -43.4;
    seeds.ROI(22).y = -71.9;
    seeds.ROI(22).z = -7.7;
    seeds.ROI(22).idx = 122;
    
    seeds.ROI(23).label = 'l-MT-2';
    seeds.ROI(23).x = -46.0;
    seeds.ROI(23).y = -66.0;
    seeds.ROI(23).z = -15.0;
    seeds.ROI(23).idx = 123;
    
    seeds.ROI(24).label = 'l-MT-3';
    seeds.ROI(24).x = -47.0;
    seeds.ROI(24).y = -68.0;
    seeds.ROI(24).z = -14.0;
    seeds.ROI(24).idx = 124;
    
    seeds.ROI(25).label = 'l-SMA';
    seeds.ROI(25).x = -5.0;
    seeds.ROI(25).y = -5.0;
    seeds.ROI(25).z = 59.0;
    seeds.ROI(25).idx = 125;
    
    seeds.ROI(26).label = 'r-IFG';
    seeds.ROI(26).x = 54.0;
    seeds.ROI(26).y = 6.0;
    seeds.ROI(26).z = 25.0;
    seeds.ROI(26).idx = 126;
    
    seeds.ROI(27).label = 'r-MFG';
    seeds.ROI(27).x = 38.0;
    seeds.ROI(27).y = 38.0;
    seeds.ROI(27).z = 20.0;
    seeds.ROI(27).idx = 127;
    
    seeds.ROI(28).label = 'r-I';
    seeds.ROI(28).x = 31.0;
    seeds.ROI(28).y = 19.0;
    seeds.ROI(28).z = 7.0;
    seeds.ROI(28).idx = 128;
   
elseif strcmp(network_label,'VAN') %%% 200
    
    seeds.network_label = 'VAN';
    
    seeds.ROI(1).label = 'r-I';
    seeds.ROI(1).x = 29.0;
    seeds.ROI(1).y = 15.0;
    seeds.ROI(1).z = -2.0;
    seeds.ROI(1).idx = 201;
    
    seeds.ROI(2).label = 'l-I';
    seeds.ROI(2).x = -43.0;
    seeds.ROI(2).y = 15.0;
    seeds.ROI(2).z = -2.0;
    seeds.ROI(2).idx = 202;
    
    seeds.ROI(3).label = 'r-TPJ';
    seeds.ROI(3).x = 57.0;
    seeds.ROI(3).y = -46.0;
    seeds.ROI(3).z = 35.0;
    seeds.ROI(3).idx = 203;
    
    seeds.ROI(4).label = 'r-SFG';
    seeds.ROI(4).x = 30.0;
    seeds.ROI(4).y = 44.0;
    seeds.ROI(4).z = 31.0;
    seeds.ROI(4).idx = 204;
    
    seeds.ROI(5).label = 'r-MidFG';
    seeds.ROI(5).x = 48.0;
    seeds.ROI(5).y = 13.0;
    seeds.ROI(5).z = 32.0;
    seeds.ROI(5).idx = 205;
    
    seeds.ROI(6).label = 'r-IFG';
    seeds.ROI(6).x = 42.0;
    seeds.ROI(6).y = 43.0;
    seeds.ROI(6).z = 2.0;
    seeds.ROI(6).idx = 206;
    
    seeds.ROI(7).label = 'r-MidTG-1';
    seeds.ROI(7).x = 55.0;
    seeds.ROI(7).y = -31.0;
    seeds.ROI(7).z = -14.0;
    seeds.ROI(7).idx = 207;
    
    seeds.ROI(8).label = 'r-MidTG-2';
    seeds.ROI(8).x = 43.0;
    seeds.ROI(8).y = -11.0;
    seeds.ROI(8).z = -14.0;
    seeds.ROI(8).idx = 208;
    
    seeds.ROI(9).label = 'l-SFG';
    seeds.ROI(9).x = -62.0;
    seeds.ROI(9).y = -50.0;
    seeds.ROI(9).z = 17.0;
    seeds.ROI(9).idx = 209;
    
    seeds.ROI(10).label = 'r-P';
    seeds.ROI(10).x = 2.0;
    seeds.ROI(10).y = -57.0;
    seeds.ROI(10).z = 57.0;
    seeds.ROI(10).idx = 210;
    
    seeds.ROI(11).label = 'r-MedFG';
    seeds.ROI(11).x = 8.0;
    seeds.ROI(11).y = -2.0;
    seeds.ROI(11).z = 69.0;
    seeds.ROI(11).idx = 211;
    
    seeds.ROI(12).label = 'r-SC';
    seeds.ROI(12).x = 6.0;
    seeds.ROI(12).y = -31.0;
    seeds.ROI(12).z = 47.0;
    seeds.ROI(12).idx = 212;
    
elseif strcmp(network_label,'SMN') %%% 300
    
    seeds.network_label = 'SMN';
    
    seeds.ROI(1).label = 'r-Cer';
    seeds.ROI(1).x = 15.0;
    seeds.ROI(1).y = -49.0;
    seeds.ROI(1).z = -29.0;
    seeds.ROI(1).idx = 301;
    
    seeds.ROI(2).label = 'l-Cer';
    seeds.ROI(2).x = -19.0;
    seeds.ROI(2).y = -50.0;
    seeds.ROI(2).z = -33.0;
    seeds.ROI(2).idx = 302;
    
    seeds.ROI(3).label = 'r-CS';
    seeds.ROI(3).x = 31.0;
    seeds.ROI(3).y = -35.0;
    seeds.ROI(3).z = 61.0;
    seeds.ROI(3).idx = 303;
    
    seeds.ROI(4).label = 'l-CS';
    seeds.ROI(4).x = -32.0;
    seeds.ROI(4).y = -35.0;
    seeds.ROI(4).z = 57.0;
    seeds.ROI(4).idx = 304;
    
    seeds.ROI(5).label = 'r-Put';
    seeds.ROI(5).x = 31.0;
    seeds.ROI(5).y = -15.0;
    seeds.ROI(5).z = 6.0;
    seeds.ROI(5).idx = 305;
    
    seeds.ROI(6).label = 'l-Put';
    seeds.ROI(6).x = -30.0;
    seeds.ROI(6).y = -16.0;
    seeds.ROI(6).z = 8.0;
    seeds.ROI(6).idx = 306;
    
    seeds.ROI(7).label = 'r-SII';
    seeds.ROI(7).x = 37.0;
    seeds.ROI(7).y = -22.0;
    seeds.ROI(7).z = 19.0;
    seeds.ROI(7).idx = 307;
    
    seeds.ROI(8).label = 'l-SII';
    seeds.ROI(8).x = -41.0;
    seeds.ROI(8).y = -25.0;
    seeds.ROI(8).z = 17.0;
    seeds.ROI(8).idx = 308;
    
    seeds.ROI(9).label = 'r-SMA';
    seeds.ROI(9).x = 4.0;
    seeds.ROI(9).y = -15.0;
    seeds.ROI(9).z = 54.0;
    seeds.ROI(9).idx = 309;
    
    seeds.ROI(10).label = 'l-SMA';
    seeds.ROI(10).x = -3.0;
    seeds.ROI(10).y = -18.0;
    seeds.ROI(10).z = 56.0;
    seeds.ROI(10).idx = 310;
    
    seeds.ROI(11).label = 'r-Tha';
    seeds.ROI(11).x = 15.0;
    seeds.ROI(11).y = -25.0;
    seeds.ROI(11).z = 2.0;
    seeds.ROI(11).idx = 311;
    
    seeds.ROI(12).label = 'l-Tha';
    seeds.ROI(12).x = -14.0;
    seeds.ROI(12).y = -24.0;
    seeds.ROI(12).z = -1.0;
    seeds.ROI(12).idx = 312;
    
    seeds.ROI(13).label = 'l-tTG';
    seeds.ROI(13).x = -42.0;
    seeds.ROI(13).y = -26.0;
    seeds.ROI(13).z = 8.0;
    seeds.ROI(13).idx = 313;
    
    seeds.ROI(14).label = 'r-sTG';
    seeds.ROI(14).x = 63.0;
    seeds.ROI(14).y = -20.0;
    seeds.ROI(14).z = 7.0;
    seeds.ROI(14).idx = 314;
    
    seeds.ROI(15).label = 'l-sTG';
    seeds.ROI(15).x = -52.0;
    seeds.ROI(15).y = 6.0;
    seeds.ROI(15).z = -11.0;
    seeds.ROI(15).idx = 315;
    
    seeds.ROI(16).label = 'r-mFG';
    seeds.ROI(16).x = 35.0;
    seeds.ROI(16).y = 39.0;
    seeds.ROI(16).z = 11.0;
    seeds.ROI(16).idx = 316;
    
    seeds.ROI(17).label = 'r-aI';
    seeds.ROI(17).x = 38.0;
    seeds.ROI(17).y = 21.0;
    seeds.ROI(17).z = 1.0;
    seeds.ROI(17).idx = 317;
    
    seeds.ROI(18).label = 'r-iPL';
    seeds.ROI(18).x = 47.0;
    seeds.ROI(18).y = -50.0;
    seeds.ROI(18).z = 50.0;
    seeds.ROI(18).idx = 318;
    
elseif strcmp(network_label,'VIS') %%% 400
    
    seeds.network_label = 'VIS';
    
    seeds.ROI(1).label = 'r-V1';
    seeds.ROI(1).x = 10.7;
    seeds.ROI(1).y = -91.8;
    seeds.ROI(1).z = 6.1;
    seeds.ROI(1).idx = 401;
    
    seeds.ROI(2).label = 'l-V1';
    seeds.ROI(2).x = -2.8;
    seeds.ROI(2).y = -100.7;
    seeds.ROI(2).z = -0.5;
    seeds.ROI(2).idx = 402;
    
    seeds.ROI(3).label = 'r-V2';
    seeds.ROI(3).x = 14.3;
    seeds.ROI(3).y = -95.7;
    seeds.ROI(3).z = 13.6;
    seeds.ROI(3).idx = 403;
    
    seeds.ROI(4).label = 'l-V2';
    seeds.ROI(4).x = -7.2;
    seeds.ROI(4).y = -98.8;
    seeds.ROI(4).z = 6.8;
    seeds.ROI(4).idx = 404;
    
    seeds.ROI(5).label = 'r-V3';
    seeds.ROI(5).x = 27.5;
    seeds.ROI(5).y = -88.5;
    seeds.ROI(5).z = 14.5;
    seeds.ROI(5).idx = 405;
    
    seeds.ROI(6).label = 'l-V3';
    seeds.ROI(6).x = -16.1;
    seeds.ROI(6).y = -92.7;
    seeds.ROI(6).z = 17.8;
    seeds.ROI(6).idx = 406;
    
    seeds.ROI(7).label = 'r-V4';
    seeds.ROI(7).x = 26.7;
    seeds.ROI(7).y = -71.4;
    seeds.ROI(7).z = -13.9;
    seeds.ROI(7).idx = 407;
    
    seeds.ROI(8).label = 'l-V4';
    seeds.ROI(8).x = -30.9;
    seeds.ROI(8).y = -76.5;
    seeds.ROI(8).z = -17.1;
    seeds.ROI(8).idx = 408;
    
    seeds.ROI(9).label = 'r-V7';
    seeds.ROI(9).x = 32.2;
    seeds.ROI(9).y = -78.4;
    seeds.ROI(9).z = 25.1;
    seeds.ROI(9).idx = 409;
    
    seeds.ROI(10).label = 'l-V7';
    seeds.ROI(10).x = -23.1;
    seeds.ROI(10).y = -78.1;
    seeds.ROI(10).z = 26.1;
    seeds.ROI(10).idx = 410;
    
    seeds.ROI(11).label = 'V2d-V3';
    seeds.ROI(11).x = 17.0;
    seeds.ROI(11).y = -101.0;
    seeds.ROI(11).z = 14.0;
    seeds.ROI(11).idx = 411;
    
    seeds.ROI(12).label = 'LO';
    seeds.ROI(12).x = 37.0;
    seeds.ROI(12).y = -90.0;
    seeds.ROI(12).z = 3.0;
    seeds.ROI(12).idx = 412;
   
elseif strcmp(network_label,'FPC') %%% 500
    
    seeds.network_label = 'FPC';
    
    seeds.ROI(1).label = 'r-AI';
    seeds.ROI(1).x = 36.0;
    seeds.ROI(1).y = 20.0;
    seeds.ROI(1).z = -2.0;
    seeds.ROI(1).idx = 501;
    
    seeds.ROI(2).label = 'l-AI';
    seeds.ROI(2).x = -37.0;
    seeds.ROI(2).y = 17.0;
    seeds.ROI(2).z = -2.0;
    seeds.ROI(2).idx = 502;
    
    seeds.ROI(3).label = 'r-DACC';
    seeds.ROI(3).x = 4.0;
    seeds.ROI(3).y = 12.0;
    seeds.ROI(3).z = 42.0;
    seeds.ROI(3).idx = 503;
    
    seeds.ROI(4).label = 'r-IPS';
    seeds.ROI(4).x = 29.0;
    seeds.ROI(4).y = -65.0;
    seeds.ROI(4).z = 42.0;
    seeds.ROI(4).idx = 504;
    
    seeds.ROI(5).label = 'l-IPS';
    seeds.ROI(5).x = -31.0;
    seeds.ROI(5).y = -62.0;
    seeds.ROI(5).z = 45.0;
    seeds.ROI(5).idx = 505;
    
    seeds.ROI(6).label = 'r-dlPFC';
    seeds.ROI(6).x = 44.0;
    seeds.ROI(6).y = 21.0;
    seeds.ROI(6).z = 34.0;
    seeds.ROI(6).idx = 506;
    
    seeds.ROI(7).label = 'l-dlPFC';
    seeds.ROI(7).x = -44.0;
    seeds.ROI(7).y = 22.0;
    seeds.ROI(7).z = 36.0;
    seeds.ROI(7).idx = 507;
    
    seeds.ROI(8).label = 'r-al/fO';
    seeds.ROI(8).x = 38.0;
    seeds.ROI(8).y = 19.0;
    seeds.ROI(8).z = 0.0;
    seeds.ROI(8).idx = 508;
    
    seeds.ROI(9).label = 'l-al/fO';
    seeds.ROI(9).x = -36.0;
    seeds.ROI(9).y = 17.0;
    seeds.ROI(9).z = 3.0;
    seeds.ROI(9).idx = 509;
   
    seeds.ROI(10).label = 'dACC/msFC';
    seeds.ROI(10).x = -2.0;
    seeds.ROI(10).y = 7.0;
    seeds.ROI(10).z = 50.0;
    seeds.ROI(10).idx = 510;
    
    seeds.ROI(11).label = 'r-aPFC';
    seeds.ROI(11).x = 28.0;
    seeds.ROI(11).y = 51.0;
    seeds.ROI(11).z = 25.0;
    seeds.ROI(11).idx = 511;
    
    seeds.ROI(12).label = 'l-aPFC';
    seeds.ROI(12).x = -28.0;
    seeds.ROI(12).y = 53.0;
    seeds.ROI(12).z = 16.0;
    seeds.ROI(12).idx = 512;
    
    seeds.ROI(13).label = 'r-TPJ';
    seeds.ROI(13).x = 54.0;
    seeds.ROI(13).y = -47.0;
    seeds.ROI(13).z = 15.0;
    seeds.ROI(13).idx = 513;
    
    seeds.ROI(14).label = 'r-midtemporal';
    seeds.ROI(14).x = 54.0;
    seeds.ROI(14).y = -32.0;
    seeds.ROI(14).z = -9.0;
    seeds.ROI(14).idx = 514;
    
    seeds.ROI(15).label = 'l-midtemporal';
    seeds.ROI(15).x = -57.0;
    seeds.ROI(15).y = -30.0;
    seeds.ROI(15).z = -11.0;
    seeds.ROI(15).idx = 515;
    
elseif strcmp(network_label,'LAN') %%% 600
    
    seeds.network_label = 'LAN';
    
    seeds.ROI(1).label = 'l-IFG';
    seeds.ROI(1).x = -48.0;
    seeds.ROI(1).y = 31.0;
    seeds.ROI(1).z = -1.0;
    seeds.ROI(1).idx = 601;
    
    seeds.ROI(2).label = 'l-MFG';
    seeds.ROI(2).x = -45.0;
    seeds.ROI(2).y = 13.0;
    seeds.ROI(2).z = 24.0;
    seeds.ROI(2).idx = 602;
    
    seeds.ROI(3).label = 'l-STS';
    seeds.ROI(3).x = -50.0;
    seeds.ROI(3).y = -54.0;
    seeds.ROI(3).z = 22.0;
    seeds.ROI(3).idx = 603;
    
    seeds.ROI(4).label = 'l-aSTG';
    seeds.ROI(4).x = -56.0;
    seeds.ROI(4).y = -12.0;
    seeds.ROI(4).z = -3.0;
    seeds.ROI(4).idx = 604;
    
    seeds.ROI(5).label = 'l-pSTG';
    seeds.ROI(5).x = -55.0;
    seeds.ROI(5).y = -48.0;
    seeds.ROI(5).z = 15.0;
    seeds.ROI(5).idx = 605;
    
    seeds.ROI(6).label = 'l-PCC-1';
    seeds.ROI(6).x = -3.0;
    seeds.ROI(6).y = -65.0;
    seeds.ROI(6).z = 29.0;
    seeds.ROI(6).idx = 606;
    
    seeds.ROI(7).label = 'l-PCC-2';
    seeds.ROI(7).x = -10.0;
    seeds.ROI(7).y = -52.0;
    seeds.ROI(7).z = 32.0;
    seeds.ROI(7).idx = 607;
    
    seeds.ROI(8).label = 'l-AG';
    seeds.ROI(8).x = -43.0;
    seeds.ROI(8).y = -70.0;
    seeds.ROI(8).z = 32.0;
    seeds.ROI(8).idx = 608;
    
    seeds.ROI(9).label = 'l-AG-latIPS';
    seeds.ROI(9).x = -34.0;
    seeds.ROI(9).y = -74.0;
    seeds.ROI(9).z = 45.0;
    seeds.ROI(9).idx = 609;
    
    seeds.ROI(10).label = 'l-aIPS-1';
    seeds.ROI(10).x = -38.0;
    seeds.ROI(10).y = -57.0;
    seeds.ROI(10).z = 37.0;
    seeds.ROI(10).idx = 610;
    
    seeds.ROI(11).label = 'l-aIPS-2';
    seeds.ROI(11).x = -46.0;
    seeds.ROI(11).y = -55.0;
    seeds.ROI(11).z = 47.0;
    seeds.ROI(11).idx = 611;
    
    seeds.ROI(12).label = 'l-PreCu';
    seeds.ROI(12).x = -2.0;
    seeds.ROI(12).y = -81.0;
    seeds.ROI(12).z = 52.0;
    seeds.ROI(12).idx = 612;
    
    seeds.ROI(13).label = 'r-AG-1';
    seeds.ROI(13).x = 47.0;
    seeds.ROI(13).y = -62.0;
    seeds.ROI(13).z = 25.0;
    seeds.ROI(13).idx = 613;
    
    seeds.ROI(14).label = 'r-AG-2';
    seeds.ROI(14).x = 43.0;
    seeds.ROI(14).y = -66.0;
    seeds.ROI(14).z = 43.0;
    seeds.ROI(14).idx = 614;
    
    seeds.ROI(15).label = 'r-AG-3';
    seeds.ROI(15).x = 59.0;
    seeds.ROI(15).y = -57.0;
    seeds.ROI(15).z = 28.0;
    seeds.ROI(15).idx = 615;
    
    seeds.ROI(16).label = 'r-latIPS';
    seeds.ROI(16).x = 34.0;
    seeds.ROI(16).y = -74.0;
    seeds.ROI(16).z = 53.0;
    seeds.ROI(16).idx = 616;
    
    seeds.ROI(17).label = 'r-PCC';
    seeds.ROI(17).x = 11.0;
    seeds.ROI(17).y = -66.0;
    seeds.ROI(17).z = 25.0;
    seeds.ROI(17).idx = 617;
    
elseif strcmp(network_label,'AUD') %%% 700
    
    seeds.network_label = 'AUD';
    
    seeds.ROI(1).label = 'r-aT';
    seeds.ROI(1).x = 45.0;
    seeds.ROI(1).y = -9.0;
    seeds.ROI(1).z = 19.0;
    seeds.ROI(1).idx = 701;
    
    seeds.ROI(2).label = 'l-aT';
    seeds.ROI(2).x = -61.0;
    seeds.ROI(2).y = -15.0;
    seeds.ROI(2).z = -3.0;
    seeds.ROI(2).idx = 702;
    
    seeds.ROI(3).label = 'r-pT';
    seeds.ROI(3).x = 45.0;
    seeds.ROI(3).y = -33.0;
    seeds.ROI(3).z = 14.0;
    seeds.ROI(3).idx = 703;
    
    seeds.ROI(4).label = 'l-pT';
    seeds.ROI(4).x = -40.0;
    seeds.ROI(4).y = -34.0;
    seeds.ROI(4).z = 10.0;
    seeds.ROI(4).idx = 704;
    
    seeds.ROI(5).label = 'r-T';
    seeds.ROI(5).x = 62.0;
    seeds.ROI(5).y = -10.0;
    seeds.ROI(5).z = -3.0;
    seeds.ROI(5).idx = 705;
    
    seeds.ROI(6).label = 'l-T';
    seeds.ROI(6).x = -60.0;
    seeds.ROI(6).y = -1.0;
    seeds.ROI(6).z = 20.0;
    seeds.ROI(6).idx = 706;
    
    seeds.ROI(7).label = 'r-MTG';
    seeds.ROI(7).x = 60.0;
    seeds.ROI(7).y = -22.0;
    seeds.ROI(7).z = 6.0;
    seeds.ROI(7).idx = 707;
    
    seeds.ROI(8).label = 'l-MTG';
    seeds.ROI(8).x = -41.0;
    seeds.ROI(8).y = -28.0;
    seeds.ROI(8).z = 6.0;
    seeds.ROI(8).idx = 708;
    
    seeds.ROI(9).label = 'r-STG';
    seeds.ROI(9).x = 54.0;
    seeds.ROI(9).y = -43.0;
    seeds.ROI(9).z = 12.0;
    seeds.ROI(9).idx = 709;
    
    seeds.ROI(10).label = 'l-STG';
    seeds.ROI(10).x = -54.0;
    seeds.ROI(10).y = -40.0;
    seeds.ROI(10).z = 10.0;
    seeds.ROI(10).idx = 710;
    
elseif strcmp(network_label,'DMN') %%% 800
    
    seeds.network_label = 'DMN';
    
    seeds.ROI(1).label = 'r-mPFC-1';
    seeds.ROI(1).x = 2.0;
    seeds.ROI(1).y = 52.6;
    seeds.ROI(1).z = 23.5;
    seeds.ROI(1).idx = 801;
    
    seeds.ROI(2).label = 'l-mPFC-1';
    seeds.ROI(2).x = -2.0;
    seeds.ROI(2).y = 50.5;
    seeds.ROI(2).z = 1.7;
    seeds.ROI(2).idx = 802;
    
    seeds.ROI(3).label = 'l-mPFC-2';
    seeds.ROI(3).x = -13.1;
    seeds.ROI(3).y = 51.5;
    seeds.ROI(3).z = 23.4;
    seeds.ROI(3).idx = 803;
    
    seeds.ROI(4).label = 'l-mPFC-3';
    seeds.ROI(4).x = 0.0;
    seeds.ROI(4).y = 59.0;
    seeds.ROI(4).z = 4.0;
    seeds.ROI(4).idx = 804;
    
    seeds.ROI(5).label = 'l-mPFC-4';
    seeds.ROI(5).x = -15.0;
    seeds.ROI(5).y = 53.0;
    seeds.ROI(5).z = 36.0;
    seeds.ROI(5).idx = 805;
    
    seeds.ROI(6).label = 'l-SFS';
    seeds.ROI(6).x = -22.2;
    seeds.ROI(6).y = 19.3;
    seeds.ROI(6).z = 51.1;
    seeds.ROI(6).idx = 806;
    
    seeds.ROI(7).label = 'l-PCC-1';
    seeds.ROI(7).x = -3.0;
    seeds.ROI(7).y = -54.0;
    seeds.ROI(7).z = 31.0;
    seeds.ROI(7).idx = 807;
    
    seeds.ROI(8).label = 'l-PCC-2';
    seeds.ROI(8).x = -5.0;
    seeds.ROI(8).y = -49.0;
    seeds.ROI(8).z = 33.0;
    seeds.ROI(8).idx = 808;
    
    seeds.ROI(9).label = 'l-RS';
    seeds.ROI(9).x = -2.0;
    seeds.ROI(9).y = -55.2;
    seeds.ROI(9).z = 10.1;
    seeds.ROI(9).idx = 809;
    
    seeds.ROI(10).label = 'r-AG-1';
    seeds.ROI(10).x = 51.0;
    seeds.ROI(10).y = -64.0;
    seeds.ROI(10).z = 32.0;
    seeds.ROI(10).idx = 810;
    
    seeds.ROI(11).label = 'r-AG-2';
    seeds.ROI(11).x = 52.0;
    seeds.ROI(11).y = -64.0;
    seeds.ROI(11).z = 27.0;
    seeds.ROI(11).idx = 811;
    
    seeds.ROI(12).label = 'r-AG-3';
    seeds.ROI(12).x = 44.0;
    seeds.ROI(12).y = -68.0;
    seeds.ROI(12).z = 53.0;
    seeds.ROI(12).idx = 812;
    
    seeds.ROI(13).label = 'l-AG-1';
    seeds.ROI(13).x = -43.0;
    seeds.ROI(13).y = -76.0;
    seeds.ROI(13).z = 35.0;
    seeds.ROI(13).idx = 813;
    
    seeds.ROI(14).label = 'l-AG-2';
    seeds.ROI(14).x = -37.0;
    seeds.ROI(14).y = -82.0;
    seeds.ROI(14).z = 35.0;
    seeds.ROI(14).idx = 814;
    
    seeds.ROI(15).label = 'l-AG-3';
    seeds.ROI(15).x = -47.0;
    seeds.ROI(15).y = -66.0;
    seeds.ROI(15).z = 45.0;
    seeds.ROI(15).idx = 815;
    
    seeds.ROI(16).label = 'l-ITG';
    seeds.ROI(16).x = -56.6;
    seeds.ROI(16).y = -25.1;
    seeds.ROI(16).z = -16.9;
    seeds.ROI(16).idx = 816;
    
    seeds.ROI(17).label = 'l-Hip';
    seeds.ROI(17).x = -23.2;
    seeds.ROI(17).y = -23.0;
    seeds.ROI(17).z = -18.0;
    seeds.ROI(17).idx = 817;
    
    seeds.ROI(18).label = 'l-mPostCentralS';
    seeds.ROI(18).x = -18.0;
    seeds.ROI(18).y = -43.0;
    seeds.ROI(18).z = 48.0;
    seeds.ROI(18).idx = 818;
    
    seeds.ROI(19).label = 'l-PostCentralS';
    seeds.ROI(19).x = -27.0;
    seeds.ROI(19).y = -42.0;
    seeds.ROI(19).z = 68.0;
    seeds.ROI(19).idx = 819;
    
    seeds.ROI(20).label = 'r-mPostCentralS';
    seeds.ROI(20).x = 11.0;
    seeds.ROI(20).y = -52.0;
    seeds.ROI(20).z = 62.0;
    seeds.ROI(20).idx = 820;
    
    seeds.ROI(21).label = 'r-mPFC-2';
    seeds.ROI(21).x = 2.0;
    seeds.ROI(21).y = 53.0;
    seeds.ROI(21).z = 13.0;
    seeds.ROI(21).idx = 821;
    
    seeds.ROI(22).label = 'l-PCC-3';
    seeds.ROI(22).x = -3.0;
    seeds.ROI(22).y = -47.0;
    seeds.ROI(22).z = 29.0;
    seeds.ROI(22).idx = 822;
    
    seeds.ROI(23).label = 'l-PCC-4';
    seeds.ROI(23).x = -44.0;
    seeds.ROI(23).y = -66.0;
    seeds.ROI(23).z = 35.0;
    seeds.ROI(23).idx = 823;
    
    seeds.ROI(24).label = 'l-PCC-5';
    seeds.ROI(24).x = -11.0;
    seeds.ROI(24).y = -52.0;
    seeds.ROI(24).z = 38.0;
    seeds.ROI(24).idx = 824;
    
    seeds.ROI(25).label = 'r-AG-4';
    seeds.ROI(25).x = 43.0;
    seeds.ROI(25).y = -66.0;
    seeds.ROI(25).z = 43.0;
    seeds.ROI(25).idx = 825;
    
    seeds.ROI(26).label = 'l-AG-4';
    seeds.ROI(26).x = -33.0;
    seeds.ROI(26).y = -80.0;
    seeds.ROI(26).z = 38.0;
    seeds.ROI(26).idx = 826;
    
    seeds.ROI(27).label = 'l-AG-5';
    seeds.ROI(27).x = -46.0;
    seeds.ROI(27).y = -71.0;
    seeds.ROI(27).z = 27.0;
    seeds.ROI(27).idx = 827;
    
end

nROI = length(seeds.ROI);

for iROI=1:nROI
   
    seeds.ROI(iROI).label = strcat(seeds.network_label,'-',seeds.ROI(iROI).label);
    
end

end