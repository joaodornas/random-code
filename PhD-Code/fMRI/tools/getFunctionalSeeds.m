function seeds = getFunctionalSeeds(network_label)

%%% coordinates in mm

if strcmp(network_label,'VAN')
   
    seeds.network_label = 'VAN';
    
    seeds.ROI(1).label = 'r-TPJ';
    seeds.ROI(1).x = 53;
    seeds.ROI(1).y = -50;
    seeds.ROI(1).z = 18;
    seeds.ROI(1).idx = 201;
    
    seeds.ROI(2).label = 'r-VFC';
    seeds.ROI(2).x = 39;
    seeds.ROI(2).y = 21;
    seeds.ROI(2).z = -3;
    seeds.ROI(2).idx = 202;
    
    seeds.ROI(3).label = 'l-TPJ';
    seeds.ROI(3).x = -53;
    seeds.ROI(3).y = -50;
    seeds.ROI(3).z = 18;
    seeds.ROI(3).idx = 203;
    
    seeds.ROI(4).label = 'l-VFC';
    seeds.ROI(4).x = -39;
    seeds.ROI(4).y = 21;
    seeds.ROI(4).z = -3;
    seeds.ROI(4).idx = 204;
    
elseif strcmp(network_label,'DAN')
    
    seeds.network_label = 'DAN';
    
    seeds.ROI(1).label = 'r-IPS';
    seeds.ROI(1).x = 26;
    seeds.ROI(1).y = -62;
    seeds.ROI(1).z = 55;
    seeds.ROI(1).idx = 101;
    
    seeds.ROI(2).label = 'r-FEF';
    seeds.ROI(2).x = 23;
    seeds.ROI(2).y = -17;
    seeds.ROI(2).z = 56;
    seeds.ROI(2).idx = 102;
    
    seeds.ROI(3).label = 'l-IPS';
    seeds.ROI(3).x = -26;
    seeds.ROI(3).y = -62;
    seeds.ROI(3).z = 55;
    seeds.ROI(3).idx = 103;
    
    seeds.ROI(4).label = 'l-FEF';
    seeds.ROI(4).x = -23;
    seeds.ROI(4).y = -17;
    seeds.ROI(4).z = 56;
    seeds.ROI(4).idx = 104;
    
elseif strcmp(network_label,'VIS')
    
    seeds.network_label = 'VIS';
    
    seeds.ROI(1).label = 'r-V1';
    seeds.ROI(1).x = 10.7;
    seeds.ROI(1).y = -91.8;
    seeds.ROI(1).z = 6.1;
    seeds.ROI(1).idx = 401;
    
    seeds.ROI(2).label = 'r-V2';
    seeds.ROI(2).x = 14.3;
    seeds.ROI(2).y = -95.7;
    seeds.ROI(2).z = 13.6;
    seeds.ROI(2).idx = 402;
    
    seeds.ROI(3).label = 'r-V3';
    seeds.ROI(3).x = 27.5;
    seeds.ROI(3).y = -88.5;
    seeds.ROI(3).z = 14.5;
    seeds.ROI(3).idx = 403;
    
    seeds.ROI(4).label = 'r-V4';
    seeds.ROI(4).x = 26.7;
    seeds.ROI(4).y = -71.4;
    seeds.ROI(4).z = -13.9;
    seeds.ROI(4).idx = 404;
    
    seeds.ROI(5).label = 'r-V7';
    seeds.ROI(5).x = 32.2;
    seeds.ROI(5).y = -78.4;
    seeds.ROI(5).z = 25.1;
    seeds.ROI(5).idx = 405;
    
    seeds.ROI(6).label = 'l-V1';
    seeds.ROI(6).x = -2.8;
    seeds.ROI(6).y = -100.7;
    seeds.ROI(6).z = -0.5;
    seeds.ROI(6).idx = 406;
    
    seeds.ROI(7).label = 'l-V2';
    seeds.ROI(7).x = -7.2;
    seeds.ROI(7).y = -98.8;
    seeds.ROI(7).z = 6.8;
    seeds.ROI(7).idx = 407;
    
    seeds.ROI(8).label = 'l-V3';
    seeds.ROI(8).x = -16.1;
    seeds.ROI(8).y = -92.7;
    seeds.ROI(8).z = 17.8;
    seeds.ROI(8).idx = 408;
    
    seeds.ROI(9).label = 'l-V4';
    seeds.ROI(9).x = -30.9;
    seeds.ROI(9).y = -76.5;
    seeds.ROI(9).z = -17.1;
    seeds.ROI(9).idx = 409;
    
    seeds.ROI(10).label = 'l-V7';
    seeds.ROI(10).x = -23.1;
    seeds.ROI(10).y = -78.1;
    seeds.ROI(10).z = 26.1;
    seeds.ROI(10).idx = 410;
    
elseif strcmp(network_label,'LAN')
    
    seeds.network_label = 'LAN';
    
    seeds.ROI(1).label = 'l-IFG';
    seeds.ROI(1).x = -48;
    seeds.ROI(1).y = 31;
    seeds.ROI(1).z = -1;
    seeds.ROI(1).idx = 601;
    
    seeds.ROI(2).label = 'l-MFG';
    seeds.ROI(2).x = -45;
    seeds.ROI(2).y = 13;
    seeds.ROI(2).z = 24;
    seeds.ROI(2).idx = 602;
    
    seeds.ROI(3).label = 'l-STS';
    seeds.ROI(3).x = -50;
    seeds.ROI(3).y = -54;
    seeds.ROI(3).z = 22;
    seeds.ROI(3).idx = 603;
    
    seeds.ROI(4).label = 'l-aSTG';
    seeds.ROI(4).x = -56;
    seeds.ROI(4).y = -12;
    seeds.ROI(4).z = -3;
    seeds.ROI(4).idx = 604;
    
    seeds.ROI(5).label = 'l-pSTG';
    seeds.ROI(5).x = -55;
    seeds.ROI(5).y = -48;
    seeds.ROI(5).z = 15;
    seeds.ROI(5).idx = 605;
   
elseif strcmp(network_label,'DMN')
    
    seeds.network_label = 'DMN';
    
    seeds.ROI(1).label = 'r-mPFC';
    seeds.ROI(1).x = 2;
    seeds.ROI(1).y = 52.6;
    seeds.ROI(1).z = 23.5;
    seeds.ROI(1).idx = 701;
    
    seeds.ROI(2).label = 'l-mPFC-1';
    seeds.ROI(2).x = -2;
    seeds.ROI(2).y = 50.5;
    seeds.ROI(2).z = 1.7;
    seeds.ROI(2).idx = 702;
    
    seeds.ROI(3).label = 'l-mPFC-2';
    seeds.ROI(3).x = -13.1;
    seeds.ROI(3).y = 51.5;
    seeds.ROI(3).z = 23.4;
    seeds.ROI(3).idx = 703;
    
    seeds.ROI(4).label = 'l-SFS';
    seeds.ROI(4).x = -22.2;
    seeds.ROI(4).y = 19.3;
    seeds.ROI(4).z = 51.1;
    seeds.ROI(4).idx = 704;
    
    seeds.ROI(5).label = 'l-PCC';
    seeds.ROI(5).x = -3;
    seeds.ROI(5).y = -54;
    seeds.ROI(5).z = 31;
    seeds.ROI(5).idx = 705;
    
    seeds.ROI(6).label = 'l-RS';
    seeds.ROI(6).x = -2;
    seeds.ROI(6).y = -55.2;
    seeds.ROI(6).z = 10.1;
    seeds.ROI(6).idx = 706;
    
    seeds.ROI(7).label = 'r-AG';
    seeds.ROI(7).x = 51;
    seeds.ROI(7).y = -64;
    seeds.ROI(7).z = 32;
    seeds.ROI(7).idx = 707;
    
    seeds.ROI(8).label = 'l-AG';
    seeds.ROI(8).x = -43;
    seeds.ROI(8).y = -76;
    seeds.ROI(8).z = 35;
    seeds.ROI(8).idx = 708;
    
    seeds.ROI(9).label = 'l-ITG';
    seeds.ROI(9).x = -56.6;
    seeds.ROI(9).y = -25.1;
    seeds.ROI(9).z = -16.9;
    seeds.ROI(9).idx = 709;
    
    seeds.ROI(10).label = 'l-HIP';
    seeds.ROI(10).x = -23.2;
    seeds.ROI(10).y = -23;
    seeds.ROI(10).z = -18;
    seeds.ROI(10).idx = 710;
    
elseif strcmp(network_label,'AUD')
    
    seeds.network_label = 'AUD';
    
    seeds.ROI(1).label = 'r-mTG';
    seeds.ROI(1).x = 60;
    seeds.ROI(1).y = -22;
    seeds.ROI(1).z = 6;
    seeds.ROI(1).idx = 801;

    seeds.ROI(2).label = 'l-mTG';
    seeds.ROI(2).x = -41;
    seeds.ROI(2).y = -28;
    seeds.ROI(2).z = 6;
    seeds.ROI(2).idx = 802;
    
    seeds.ROI(3).label = 'r-sTG';
    seeds.ROI(3).x = 54;
    seeds.ROI(3).y = -43;
    seeds.ROI(3).z = 12;
    seeds.ROI(3).idx = 803;
    
    seeds.ROI(4).label = 'l-sTG';
    seeds.ROI(4).x = -54;
    seeds.ROI(4).y = -40;
    seeds.ROI(4).z = 10;
    seeds.ROI(4).idx = 804;

elseif strcmp(network_label,'FPC')
    
    seeds.network_label = 'FPC';
    
    seeds.ROI(1).label = 'r-AI';
    seeds.ROI(1).x = 36;
    seeds.ROI(1).y = 20;
    seeds.ROI(1).z = -2;
    seeds.ROI(1).idx = 501;
    
    seeds.ROI(2).label = 'l-AI';
    seeds.ROI(2).x = -37;
    seeds.ROI(2).y = 17;
    seeds.ROI(2).z = -2;
    seeds.ROI(2).idx = 502;
    
    seeds.ROI(3).label = 'r-dACC';
    seeds.ROI(3).x = 4;
    seeds.ROI(3).y = 12;
    seeds.ROI(3).z = 42;
    seeds.ROI(3).idx = 503;
    
elseif strcmp(network_label,'SMN')
    
    seeds.network_label = 'SMN';
    
    seeds.ROI(1).label = 'l-Cbllm';
    seeds.ROI(1).x = -19;
    seeds.ROI(1).y = -50;
    seeds.ROI(1).z = -33;
    seeds.ROI(1).idx = 301;

    seeds.ROI(2).label = 'l-CS';
    seeds.ROI(2).x = -32;
    seeds.ROI(2).y = -35;
    seeds.ROI(2).z = 57;
    seeds.ROI(2).idx = 302;
    
    seeds.ROI(3).label = 'l-Put';
    seeds.ROI(3).x = -30;
    seeds.ROI(3).y = -16;
    seeds.ROI(3).z = 8;
    seeds.ROI(3).idx = 303;
    
    seeds.ROI(4).label = 'l-SII';
    seeds.ROI(4).x = -41;
    seeds.ROI(4).y = -25;
    seeds.ROI(4).z = 17;
    seeds.ROI(4).idx = 304;
    
    seeds.ROI(5).label = 'l-SMA';
    seeds.ROI(5).x = -3;
    seeds.ROI(5).y = -18;
    seeds.ROI(5).z = 56;
    seeds.ROI(5).idx = 305;
    
    seeds.ROI(6).label = 'l-Thal';
    seeds.ROI(6).x = -14;
    seeds.ROI(6).y = -24;
    seeds.ROI(6).z = -1;
    seeds.ROI(6).idx = 306;
    
    seeds.ROI(7).label = 'r-Cbllm';
    seeds.ROI(7).x = 15;
    seeds.ROI(7).y = -49;
    seeds.ROI(7).z = -29;
    seeds.ROI(7).idx = 307;
    
    seeds.ROI(8).label = 'r-CS';
    seeds.ROI(8).x = 31;
    seeds.ROI(8).y = -35;
    seeds.ROI(8).z = 61;
    seeds.ROI(8).idx = 308;
    
    seeds.ROI(9).label = 'r-Put';
    seeds.ROI(9).x = 31;
    seeds.ROI(9).y = -15;
    seeds.ROI(9).z = 6;
    seeds.ROI(9).idx = 309;
    
    seeds.ROI(10).label = 'r-SII';
    seeds.ROI(10).x = 37;
    seeds.ROI(10).y = -22;
    seeds.ROI(10).z = 19;
    seeds.ROI(10).idx = 310;
    
    seeds.ROI(11).label = 'r-SMA';
    seeds.ROI(11).x = 4;
    seeds.ROI(11).y = -15;
    seeds.ROI(11).z = 54;
    seeds.ROI(11).idx = 311;
    
    seeds.ROI(12).label = 'r-Thal';
    seeds.ROI(12).x = 15;
    seeds.ROI(12).y = -25;
    seeds.ROI(12).z = 2;
    seeds.ROI(12).idx = 312;
    
end

end