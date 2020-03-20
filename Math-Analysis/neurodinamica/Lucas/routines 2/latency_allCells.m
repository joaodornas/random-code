%load variable declared below from workspace with file dialog box. 
%variable is a character array which contains cell names 
load ('population_data','cells','complex_i');
complex_cells=cells(complex_i);
variable='complex_cells';

%declare input parameters and file path for saving
method='maxLikelihood';
savepath='/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/latency/';
psthConvolution='convolved';
psthBinSize='1';
sigma='3';

%create matrices
ncells=eval(['size(' variable ',1)']);
latency_allCells_maxLikelihood=zeros(ncells,36);

% calculate latency for every condition and every cell, and save in a matrix
% latency_allCells_method, where rows are cells and columns are conditions.
% also calculate evoked response significance and replace latency values for
% NaN in conditions where there is no evoked response

for k=[1:4,6:12,14:ncells]

    disp(['computing cell #' num2str(k) '...'])

    cellname=char(eval([variable '(k)']));
    filename=['psth_' cellname '_' psthConvolution '_sigma3.mat'];
    load(filename);

    parameters.analysis_period = 250;
    parameters.analysis_startTime = 20;
    parameters.alpha = 0.05;
    nconditions=parameters.nconditions;

    [h]=evokedresponse(spike_times,stimIds,parameters);

    i_nonsig=find(h==0)';
    latency_allCells_maxLikelihood(k,i_nonsig)=NaN;

    for i=1:nconditions

        if ~isnan(latency_allCells_maxLikelihood(k,i))

            psth_i=eval(['psth.condition' num2str(i)]);
            l=maxLikelihood_latency(psth_i,parameters);

            if isempty(l)==1
                latency_allCells_maxLikelihood(k,i)=NaN;
            else latency_allCells_maxLikelihood(k,i)=l;
            end

        end

    end

    clear l cellname psth spike_times stimIds psth_i filename h nconditions
    save([savepath 'latency_' method '_psth_' psthConvolution '_binsize' psthBinSize '_sigma3.mat']);

end

% %save file in folder specified in savepath
% save([savepath 'latency_' method '_psth_' psthConvolution '_binsize' psthBinSize '.mat']);