

cells_data_file = get_all_cells_FB;

for iCell=1:length(cells_data_file)
%for iCell=1:1
    
    nVideos = length(cells_data_file(iCell).registro_video);
    
    for iVideo=1:nVideos
    %for iVideo=1:1
        
        registro = char(strcat(cells_data_file(iCell).registro_video(iVideo).datafile,'.mat'));
        
        name = char(strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile));
        
        name = strrep(name,'_','');
        
        mkdir(name);
        
        latency_file = strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile,'-','latency','.mat');

        latency = load(latency_file);

        Spass = load(char(registro));

        spike_times = Spass.spike_times ./ 32000;

        start_time = 500;
        if iCell == 1, end_time = 3500; else end_time = 9500; end
        total_time = end_time - start_time;
        
        nConditions = max(Spass.stimIds);

        opts.shift_cost = [0 2.^(-4:9)];
        resolution = ( 1./opts.shift_cost );
        resolution(1) = [];
        nResolutions = length(resolution);

        for iCondition=1:nConditions;
        %for iCondition=1:1
            
            disp(strcat('Condition:',num2str(iCondition)));  

            labels = find(Spass.stimIds == iCondition); 

            trials_spikes = spike_times(labels(:),:);
            
            trials_spikes = trials_spikes - latency.latency_time(iCondition);

            nTrials(iCondition) = size(trials_spikes,1);

            for iResolution=1:nResolutions
                
                disp(strcat('Resolution:',num2str(resolution(iResolution))));
                
                sigma = resolution(iResolution)*10^3; 
                
                all_trains_conv = zeros(nTrials(iCondition),total_time);
                all_trains_delta = zeros(nTrials(iCondition),total_time);

                for iTrial=1:nTrials(iCondition)
                    
                    if mod(iTrial,10) == 0, 
                        disp(strcat('Trial:',num2str(iTrial))); 
                    end                   
                    
                       [spike_train, delta_functions] = getSpikesConvolved(trials_spikes(iTrial,:),sigma,total_time,start_time,end_time);

                    if iTrial == 1

                        f = figure;
                        
                        plot(delta_functions,'b');
                        hold on
                        
                        plot(spike_train,'r');
                        
                        printPlot = true;

                    end
                    
                       all_trains_conv(iTrial,:) = spike_train;
                       all_trains_delta(iTrial,:) = delta_functions;

                end
                
                reliability(iResolution).reliability = getReliability(all_trains_conv);

                reliability(iResolution).all_trains_conv = all_trains_conv;
                
                reliability(iResolution).all_trains_delta = all_trains_delta;
                
                if printPlot
                
                    title(strcat('Reliability:',num2str(reliability(iResolution).reliability),'-','Resolution:',num2str(resolution(iResolution))));

                    print(f,'-depsc',strcat(name,'/',name,'-','Condition','-',int2str(iCondition),'-',int2str(iResolution),'-', 'Reliability-',num2str(reliability(iResolution).reliability),'-','Resolution-',num2str(resolution(iResolution)),'.eps'));
                        
                    close all
                    
                    printPlot = false;
                    
                end
                
                clear all_trains_conv
                clear all_trains_delta

            end

            Condition(iCondition).reliability = reliability;
            
            clear reliability
            
        end
        
        disp('Conditions Together');
        
        nSpikes = size(spike_times,2);
        
        all_trials = zeros(nTrials(1)+nTrials(2),nSpikes);
        
        for iCondition=1:2
            
            labels = find(Spass.stimIds == iCondition); 

            trials_spikes = spike_times(labels(:),:);
            
            trials_spikes = trials_spikes - latency.latency_time(iCondition);
            
            if iCondition == 2
                
                 nBins = (end_time - start_time);
                 
                 for iTrial=1:nTrials(iCondition)
                     
                     spikes = trials_spikes(iTrial,:);

                     spikes = inverteTrilha(spikes,start_time,nBins);
                     
                     spikes = padarray(spikes',[size(trials_spikes,2)-size(spikes,2) 0],'symmetric','pos')';
                     
                     trials_spikes_invertido(iTrial,:) = spikes;
                     
                 end
                 
            end
            
            if iCondition == 1, start_count = 1; else start_count = nTrials(1) + 1; end
            if iCondition == 1, end_count = nTrials(1); else end_count = nTrials(1) + nTrials(2); end
            
            all_trials(start_count:end_count,1:nSpikes) = trials_spikes;
            
            if iCondition == 1
                
                all_trials_invertido(start_count:end_count,1:nSpikes) = trials_spikes;
                
            else
                
                all_trials_invertido(start_count:end_count,1:nSpikes) = trials_spikes_invertido;
                
            end
            
        end
        
        clear trials_spikes
        clear trials_spikes_invertido
        
        for iResolution=1:nResolutions
                
            disp(strcat('Resolution:',num2str(resolution(iResolution))));
                
            sigma = resolution(iResolution)*10^3; 

            for iTrial=1:nTrials(1)             
                    
                [spike_train, delta_functions] = getSpikesConvolved(all_trials(iTrial,:),sigma,total_time,start_time,end_time);
                
                all_trains_conv(iTrial,:) = spike_train;
                all_trains_delta(iTrial,:) = delta_functions;
                
            end
            
            for iTrial=1:nTrials(2)

                [spike_train_invertido, delta_functions_invertido] = getSpikesConvolved(all_trials_invertido(iTrial,:),sigma,total_time,start_time,end_time);

                all_trains_conv_invertido(iTrial,:) = spike_train_invertido;
                all_trains_delta_invertido(iTrial,:) = delta_functions_invertido;
                
            end
            
            reliability_WO_inversion(iResolution).reliability = getReliability(all_trains_conv);

            reliability_WO_inversion(iResolution).all_trains_conv = all_trains_conv;
                
            reliability_WO_inversion(iResolution).all_trains_delta = all_trains_delta;
            
            reliability_WITH_inversion(iResolution).reliability = getReliability(all_trains_conv_invertido);

            reliability_WITH_inversion(iResolution).all_trains_conv = all_trains_conv_invertido;
                
            reliability_WITH_inversion(iResolution).all_trains_delta = all_trains_delta_invertido;
            
        end
        
        clear all_trials
        clear all_trials_invertido
        
        clear all_trains_conv
        clear all_trains_delta
        clear all_trains_conv_invertido
        clear all_trains_delta_invertido
        
        Together.reliability_WO_inversion = reliability_WO_inversion;
        Together.reliability_WITH_inversion = reliability_WITH_inversion;
        
        colors{1} = 'b-o';
        colors{2} = 'b--o';
        colors{3} = 'y-o';
        colors{4} = 'r--d';
        colors{5} = 'r-d';
        
        markerSize(1) = 4;
        markerSize(2) = 4;
        markerSize(3) = 4;
        markerSize(4) = 10;
        markerSize(5) = 10;
        
        lineWidth = 1;
        
        f = figure;
        
        for iCondition=1:nConditions
            
            for iResolution=1:nResolutions
                
                real(iResolution) = Condition(iCondition).reliability(iResolution).reliability;
        
            end
            
            plot(1:length(real),real,colors{iCondition},'LineWidth',lineWidth,'MarkerSize',markerSize(iCondition));
           
            hold on;
            
        end
        
        for iResolution=1:nResolutions
            
            real_WO(iResolution) = reliability_WO_inversion(iResolution).reliability;
            
            real_WITH(iResolution) = reliability_WITH_inversion(iResolution).reliability;
            
        end
        
        plot(1:length(real),real_WO,colors{4},'LineWidth',lineWidth,'MarkerSize',markerSize(4));
        
        hold on;
        
        plot(1:length(real),real_WITH,colors{5},'LineWidth',lineWidth,'MarkerSize',markerSize(5));

        ylim([0 1]);
        xlim([0 14]);
        set(gca,'XTick',0:1:14);
        set(gca,'xticklabel',[0,resolution],'FontSize',5);
        
        xlabel('Resolutions (s)');
        ylabel('Reliabilities');
        
        legend({'Forward','Backward','Spontaneous Activity','Forward and Backward (without inversion)','Forward and Backward (with inversion)'});
        
        title(name);
        
        print(f,'-djpeg',strcat(name,'-reliability.jpeg'));
        print(f,'-depsc',strcat(name,'-reliability.eps'));
        
        close all;
        
        clear reliability;
            
        
        save(strcat('ForBack-',name,'-reliability.mat'),'Condition','Together');
        
%         clear Condition
%         clear Together
        
    end
    
end
                
   





