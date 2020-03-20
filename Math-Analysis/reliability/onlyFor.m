function onlyFor(registro,j)

Spass = load(char(registro));

bandwidths = 1:500;

HSIZE = 1:500;

forwardlabels = find(Spass.stimIds == 1);

spike_times = Spass.spike_times ./ 32000;

forwardtrials = spike_times(forwardlabels(:),:);

nTrials_for = size(forwardtrials,1);

start_time = 500;

end_time = 9500;

mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/');

Filter_Dimension = 1;

name = char(Spass.cellname);

bar = '/';

filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-Forward-HSIZE-',int2str(HSIZE(j)),'-',int2str(Filter_Dimension),'D.mat');


[a, b] = getAllReliabilities(HSIZE(j),bandwidths,forwardtrials,nTrials_for,start_time,end_time,Filter_Dimension);

TV_for = a;
real_for = b;

real_for(isnan(real_for)) = 0;

max_for_real = max(real_for);

max_resolution_for = find(real_for==max_for_real);

infoFor = struct('name',name,'real_for',real_for,'TV_for',TV_for,'max_for_real',max_for_real,'max_resolution_for',max_resolution_for,'Filter_Dimension',Filter_Dimension);

save(filename, 'infoFor','-v7.3');

end

