function reliabilityOptimal(start_registro,end_registro,OS)

tic

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
    Spass(r) = load(char(registro(r)));

end

bandwidths = 1:9000;

e = 0.01;

M = bandwidths .* sqrt( -2 * log ( e ) ) ;

HSIZE = round ( ( 2.*M ) + 1 );

D = 1;
    
%getReal(Spass,start_registro,end_registro,HSIZE,bandwidths,'one_by_one',D,OS);
getReal(Spass,start_registro,end_registro,HSIZE,bandwidths,'all_together',D,OS);

%'one_by_one' means one HSIZE per time;
%'all_together' means all HSIZE at the same time;

function getReal(Spass,start_registro,end_registro,HSIZE,bandwidths,how_to_do,D,OS)
        
    nConditions = 3;

    start_time = 500;

    end_time = 9500;
    
    if strcmp(OS,'Mac')
        
        mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/',how_to_do,'/');
        
        bar = '/';
        
    elseif strcmp(OS,'Win')
        
        mainPath = strcat('Z:\DATA\Forward-Backward\reliabilityOptimal\',how_to_do,'\');
        
        bar = '\';
        
    end
    
    for h=start_registro:end_registro

        
        %%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(strcat('Read trials...',char(Spass(h).cellname)));

        name = char(Spass(h).cellname);
        
        mkdir(mainPath,name);
        
        forwardlabels = find(Spass(h).stimIds == 1); 
        
        backwardlabels = find(Spass(h).stimIds == 2); 
        
        AE_labels = find(Spass(h).stimIds == 3); 

        spike_times = Spass(h).spike_times ./ 32000;

        forwardtrials = spike_times(forwardlabels(:),:);

        backwardtrials = spike_times(backwardlabels(:),:);
        
        AE_trials = spike_times(AE_labels(:),:);
        
        forwardlabels = [];
        backwardlabels = [];
        AE_labels = [];
 
        nTrials_for = size(forwardtrials,1);
        
        nTrials_back = size(backwardtrials,1);
        
        nTrials_ae = size(AE_trials,1);
        
        
        for Filter_Dimension=1:D
        
                
            if strcmp(how_to_do,'one_by_one')

                for j=1:length(HSIZE)
                    
                    filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-Forward-HSIZE-',int2str(HSIZE(j)),'-',int2str(Filter_Dimension),'D.mat');
   
                    if exist(filename,'file') == 0 
                   
                        [a, b] = getAllReliabilities(HSIZE(j),bandwidths,forwardtrials,nTrials_for,start_time,end_time,Filter_Dimension,how_to_do);

                        TV_for = a;
                        real_for = b;

                        real_for(isnan(real_for)) = 0;

                        max_for_real = max(real_for);

                        max_resolution_for = find(real_for==max_for_real);

                        infoFor = struct('name',name,'real_for',real_for,'TV_for',TV_for,'max_for_real',max_for_real,'max_resolution_for',max_resolution_for,'Filter_Dimension',Filter_Dimension);

                        save(filename, 'infoFor','-v7.3');

                        infoFor = [];
                        real_for = [];
                        TV_for = [];
                    
                    end

                end

            elseif strcmp(how_to_do,'all_together')
                
                filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-Forward-HSIZE-min-',int2str(min(HSIZE)),'-max-',int2str(max(HSIZE)),'-',int2str(Filter_Dimension),'D.mat');

                if exist(filename,'file') == 0
                    
                    [TV_for, real_for] = getAllReliabilities(HSIZE,bandwidths,forwardtrials,nTrials_for,start_time,end_time,Filter_Dimension,how_to_do);

                    real_for(isnan(real_for)) = 0;

                    max_for_real = max(real_for);

                    max_resolution_for = find(real_for==max_for_real);

                    infoFor = struct('name',name,'real_for',real_for,'TV_for',TV_for,'max_for_real',max_for_real,'max_resolution_for',max_resolution_for,'Filter_Dimension',Filter_Dimension);

                    save(filename, 'infoFor','-v7.3');

                    infoFor = [];
                    real_for = [];
                    TV_for = [];
                
                end

            end


           if strcmp(how_to_do,'one_by_one')

             for j=1:length(HSIZE)
                 
                    filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-Backward-HSIZE-',int2str(HSIZE(j)),'-',int2str(Filter_Dimension),'D.mat');

                    if exist(filename,'file') == 0
                    
                        [a, b] = getAllReliabilities(HSIZE(j),bandwidths,backwardtrials,nTrials_back,start_time,end_time,Filter_Dimension,how_to_do);

                        TV_back = a;
                        real_back = b;

                        real_back(isnan(real_back)) = 0;

                        max_back_real = max(real_back);

                        max_resolution_back = find(real_back==max_back_real);

                        infoBack = struct('name',name,'real_back',real_back,'TV_back',TV_back,'max_back_real',max_back_real,'max_resolution_back',max_resolution_back,'Filter_Dimension',Filter_Dimension);

                        save(filename, 'infoBack','-v7.3');

                        infoBack = [];
                        real_back = [];
                        TV_back = [];
                    
                    end

                end

            elseif strcmp(how_to_do,'all_together')
                
                filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-Backward-HSIZE-min-',int2str(min(HSIZE)),'-max-',int2str(max(HSIZE)),'-',int2str(Filter_Dimension),'D.mat');

                if exist(filename,'file') == 0
                    
                    [TV_back, real_back] = getAllReliabilities(HSIZE,bandwidths,backwardtrials,nTrials_back,start_time,end_time,Filter_Dimension,how_to_do);

                    real_back(isnan(real_back)) = 0;

                    max_back_real = max(real_back);

                    max_resolution_back = find(real_back==max_back_real);

                    infoBack = struct('name',name,'real_back',real_back,'TV_back',TV_back,'max_back_real',max_back_real,'max_resolution_back',max_resolution_back,'Filter_Dimension',Filter_Dimension);

                    save(filename, 'infoBack','-v7.3');

                    infoBack = [];
                    real_back = [];
                    TV_back = [];
                
                end

            end

            

           if strcmp(how_to_do,'one_by_one')

                for j=1:length(HSIZE)

                    filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-AE-HSIZE-',int2str(HSIZE(j)),'-',int2str(Filter_Dimension),'D.mat');
                    
                    if exist(filename,'file') == 0
                        
                        [a, b] = getAllReliabilities(HSIZE(j),bandwidths,AE_trials,nTrials_ae,start_time,end_time,Filter_Dimension,how_to_do);

                        TV_AE = a;
                        real_AE = b;

                        real_AE(isnan(real_AE)) = 0;

                        max_AE_real = max(real_AE);

                        max_resolution_AE = find(real_AE==max_AE_real);

                        info_AE = struct('name',name,'real_AE',real_AE,'TV_AE',TV_AE,'max_AE_real',max_AE_real,'max_resolution_AE',max_resolution_AE,'Filter_Dimension',Filter_Dimension);

                        save(filename,'info_AE','-v7.3');

                        info_AE = [];
                        real_AE = [];
                        TV_AE = [];
                    
                    end

                end

            elseif strcmp(how_to_do,'all_together')
                
                filename = strcat(mainPath,name,bar,name,'-reliabilityOptimal-AE-HSIZE-min-',int2str(min(HSIZE)),'-max-',int2str(max(HSIZE)),'-',int2str(Filter_Dimension),'D.mat');

                if exist(filename,'file') == 0
                    
                    [TV_AE, real_AE] = getAllReliabilities(HSIZE,bandwidths,backwardtrials,nTrials_back,start_time,end_time,Filter_Dimension,how_to_do);

                    real_AE(isnan(real_AE)) = 0;

                    max_AE_real = max(real_AE);

                    max_resolution_AE = find(real_AE==max_AE_real);

                    info_AE = struct('name',name,'real_AE',real_AE,'TV_AE',TV_AE,'max_AE_real',max_AE_real,'max_resolution_AE',max_resolution_AE,'Filter_Dimension',Filter_Dimension);

                    save(filename,'info_AE','-v7.3');

                    info_AE = [];
                    real_AE = [];
                    TV_AE = [];
                
                end

           end 
     
        end
        
        forwardtrials = [];
        backwardtrials = [];
        AE_trials = [];
        
        Spass(h).spike_times = [];
        Spass(h).stimIds = [];
        spike_times = [];
        
        disp(strcat('Finish trials...',char(Spass(h).cellname)));
        
    end
   
end

toc

end

