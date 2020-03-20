

analysis = 'Functional-Parcels-Pop-Map';

%settings_jan_0805;
%settings_elena_2905;

threshold = 3;
pvalue = 0.05;
volume = 22;
smooth = 6;

prefix = strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis);

netseed{1} = 'DAN-l-FEF';
netseed{2} = 'DAN-r-FEF';
netseed{3} = 'DAN-l-IPS';
netseed{4} = 'DAN-r-IPS';
netseed{5} = 'VAN-l-TPJ';
netseed{6} = 'VAN-r-TPJ';
netseed{7} = 'VAN-l-VFC';
netseed{8} = 'VAN-r-VFC';

netseed{9} = 'LAN-l-IFG';
netseed{10} = 'LAN-l-MFG';
netseed{11} = 'LAN-l-aSTG';
netseed{12} = 'LAN-l-pSTG';
netseed{13} = 'VIS-l-V1';
netseed{14} = 'VIS-l-V2';
netseed{15} = 'VIS-l-V3';
netseed{16} = 'VIS-l-V4';
netseed{17} = 'VIS-l-V7';
netseed{18} = 'VIS-r-V1';
netseed{19} = 'VIS-r-V2';
netseed{20} = 'VIS-r-V3';
netseed{21} = 'VIS-r-V4';
netseed{22} = 'VIS-r-V7';

for iseed=1:length(netseed)
    
    inputimg = strcat(prefix,'-',netseed{iseed});
    outputimg = strcat(inputimg,'-clu');
    system (sprintf('/usr/local/fsl/bin/cluster --in=%s --zthresh=%s --othresh=%s --dlh=%s --volume=%s --pthresh=%s',inputimg,int2str(threshold),outputimg,int2str(smooth),int2str(volume),num2str(pvalue)));

end
