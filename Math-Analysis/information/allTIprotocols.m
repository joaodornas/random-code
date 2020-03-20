function allTIprotocols


%% GRATINGS - CENTER SURROUND - 44PX

protocols = 'cc-44px-2sizes-g-Protocols.txt';
start_stimulus_time = 500;
end_stimulus_time = 2500;

folder = '/Users/joaodornas/Documents/Dropbox/_Research/_DATA/Center-Surround/information/';
%folder = 'C:\Users\Dornas\Documents\Dropbox\_Research\_DATA\Center-Surround\information\';

%mutualInfo(protocols,start_stimulus_time,end_stimulus_time,folder);

%allTIPlotprotocols(protocols,start_stimulus_time,end_stimulus_time,folder);

%% NATURAL SCENES - CENTER SURROUNG - 44PX

protocols = 'cc-44px-2sizes-v-Protocols.txt';
start_stimulus_time = 500;
end_stimulus_time = 9500;

folder = '/Users/joaodornas/Documents/Dropbox/_Research/_DATA/Center-Surround/information/';
%folder = 'C:\Users\Dornas\Documents\Dropbox\_Research\_DATA\Center-Surround\information\';

%mutualInfo(protocols,start_stimulus_time,end_stimulus_time,folder);

%allTIPlotprotocols(protocols,start_stimulus_time,end_stimulus_time,folder);

%% NATURAL SCENES - FORWARD-BACKWARD 

protocols = 'FB-v-protocols.txt';
start_stimulus_time = 500;
end_stimulus_time = 9500;

folder = '/Users/joaodornas/Documents/Dropbox/_Research/_DATA/Forward-Backward/information/';
%folder = 'C:\Users\Dornas\Documents\Dropbox\_Research\_DATA\Forward-Backward\information\';

%mutualInfo(protocols,start_stimulus_time,end_stimulus_time,folder);

%allTIPlotprotocols(protocols,start_stimulus_time,end_stimulus_time,folder);

% GRATINGS - FORWARD-BACKWARD

folder = '/Users/joaodornas/Documents/Dropbox/_Research/_DATA/Forward-Backward/information/';
%folder = 'C:\Users\Dornas\Documents\Dropbox\_Research\_DATA\Forward-Backward\information\';

registros = importdata(protocols);

for p=1:length(registros)
    
    %[start_stimulus_time, end_stimulus_time, which_one] = DTC_forback(char(registros{p}));

    %if ~strcmp(which_one,'none')
        
        mutualInfo('blkdtc002a01_2b.mat',start_stimulus_time,end_stimulus_time,folder);
        
        %allTIPlotprotocols(which_one,start_stimulus_time,end_stimulus_time,folder);
        
    %end
    
end     


end

