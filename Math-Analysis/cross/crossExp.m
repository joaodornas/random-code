function crossExp(registro,qual_cross)

tic

disp('BEGIN');

nConditions = 2;

start_time = 500;

end_time = 9500;

latency_end_time = start_time + 300;

bin_size = 1;

p = 10000;

nBins = (end_time - start_time)/bin_size;

%%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read trials...');
     
    Spass = load(char(registro));
    
    for i=1:nConditions
   
        trials_labels(i).labels = find(Spass.stimIds == i); 
    
    end
    
    spike_times = Spass.spike_times./ 32000;
    
    forwardtrials = spike_times(trials_labels(1).labels(:),:);
    
    forwardtrials = reshape(forwardtrials.',[],1);
    
    forwardtrials = forwardtrials.';
    
    forwardtrials = sort(forwardtrials);
    
    forwardtrials = forwardtrials(forwardtrials>0);
    
    forwardtrials = forwardtrials(forwardtrials>=0 & forwardtrials<=(end_time/1000));
    
    %forwardtrialskernel = forwardtrials(forwardtrials>=0 & forwardtrials<=(latency_end_time/1000));
        
    backwardtrials = spike_times(trials_labels(2).labels(:),:);
        
    backwardtrials = reshape(backwardtrials.',[],1);

    backwardtrials = backwardtrials.';
    
    backwardtrials = sort(backwardtrials);
    
    backwardtrials = backwardtrials(backwardtrials>0);    
    
    backwardtrials = backwardtrials(backwardtrials>=0 & backwardtrials<=(end_time/1000));
    
    %backwardtrialskernel = backwardtrials(backwardtrials>=0 & backwardtrials<=(latency_end_time/1000));
    
    protocol = registro(2:end-4);

    latencyFor = getLatencyForBack(protocol,1);
    
    latencyBack = getLatencyForBack(protocol,2);
    
    if qual_cross == 1
        
        autoforcorr = xcorr(forwardtrials,forwardtrials,'none');
        
        autobackcorr = xcorr(backwardtrials,backwardtrials,'none');
    
        forbackcorr = xcorr(forwardtrials,backwardtrials,'none');
        
    elseif qual_cross == 2
    
        if (length(forwardtrials) > 0) && (length(backwardtrials) > 0)
            
            forbackcorr = correlogram(forwardtrials,backwardtrials);
            
        else
            
            forbackcorr = zeros(1,end_time);
            
        end
        
    elseif qual_cross == 3
        
        autoforcorr = getCorrelation(forwardtrials,forwardtrials);
        
        autobackcorr = getCorrelation(backwardtrials,backwardtrials);
    
        forbackcorr = getCorrelation(forwardtrials,backwardtrials);
        
    end
    
    a_e_for = forwardtrials(forwardtrials>=(0/1000) & forwardtrials<(start_time/1000));
    
    a_e_back = backwardtrials(backwardtrials>=(0/1000) & backwardtrials<(start_time/1000));

    after_ae_for = forwardtrials(forwardtrials>=(start_time/1000) & forwardtrials<(end_time/1000)) - latencyFor;
    
    after_ae_back = backwardtrials(backwardtrials>=(start_time/1000) & backwardtrials<(end_time/1000)) - latencyBack;
       
    forWOlatency = [];
    
    backWOlatency = [];
    
    forWOlatency = [a_e_for, after_ae_for];
    
    backWOlatency = [a_e_back, after_ae_back];
    
    forWOlatency = forWOlatency(forWOlatency>0);
    
    backWOlatency = backWOlatency(backWOlatency>0);
    
    if qual_cross == 1
            
        forbackWOlatencycorr = xcorr(forWOlatency,backWOlatency,'none');
        
    elseif qual_cross == 2
        
        if (length(forWOlatency) > 0) && (length(backWOlatency) > 0)
    
            forbackWOlatencycorr = correlogram(forWOlatency,backWOlatency);
            
        else
            
            forbackWOlatencycorr = zeros(1,end_time);
            
        end
        
    elseif qual_cross == 3
        
        forbackWOlatencycorr = getCorrelation(forWOlatency,backWOlatency);
        
    end
    
    backWOlatencyInverted = backWOlatency(backWOlatency>=(start_time/1000) & backWOlatency<=(end_time/1000));
    
    backWOlatencyInverted = invertSpikeTrain(backWOlatencyInverted,start_time,nBins,bin_size);
    
    backWOlatencyInverted = backWOlatencyInverted.spikes;
    
    backWOlatencyInverted = backWOlatencyInverted(backWOlatencyInverted>0);
    
    backInverted = backwardtrials(backwardtrials>=(start_time/1000) & backwardtrials<=(end_time/1000));
    
    backInverted = invertSpikeTrain(backInverted,start_time,nBins,bin_size);
    
    backInverted = backInverted.spikes;
    
    backInverted = backInverted(backInverted>0);
    
    
    if qual_cross == 1
        
        forinvertedbackWOlatencycorr = xcorr(forWOlatency,backWOlatencyInverted,'none');
        
        forwardinvertedbackcorr = xcorr(forwardtrials,backInverted,'none');
        
    elseif qual_cross == 2
        
        if (length(forWOlatency) > 0) && (length(backWOlatencyInverted) > 0)
        
            forinvertedbackWOlatencycorr = correlogram(forWOlatency,backWOlatencyInverted);
            
        else
            
            forinvertedbackWOlatencycorr = zeros(1,end_time);
            
        end
        
       if (length(forwardtrials) > 0) && (length(backInverted) > 0)
        
            forwardinvertedbackcorr = correlogram(forwardtrials,backInverted);
            
        else
            
            forwardinvertedbackcorr = zeros(1,end_time);
            
       end
        
    elseif qual_cross == 3
        
        forinvertedbackWOlatencycorr = getCorrelation(forWOlatency,backWOlatencyInverted);
        
        forwardinvertedbackcorr = getCorrelation(forwardtrials,backInverted);
        
    end
    
    name = registro;
    
    name = name(1:end-4);
    
    if qual_cross == 1
    
        correlation = struct('name',name,'autoforcorr',autoforcorr,'autobackcorr',autobackcorr,'forbackcorr',forbackcorr,'forbackWOlatencycorr',forbackWOlatencycorr,'forinvertedbackWOlatencycorr',forinvertedbackWOlatencycorr,'forwardinvertedcorr',forwardinvertedbackcorr,'latencyFor',latencyFor,'latencyBack',latencyBack);
    
    elseif qual_cross == 2
        
        correlation = struct('name',name,'forbackcorr',forbackcorr,'forbackWOlatencycorr',forbackWOlatencycorr,'forinvertedbackWOlatencycorr',forinvertedbackWOlatencycorr,'forwardinvertedcorr',forwardinvertedbackcorr,'latencyFor',latencyFor,'latencyBack',latencyBack);
            
    elseif qual_cross == 3
        
        correlation = struct('name',name,'autoforcorr',autoforcorr,'autobackcorr',autobackcorr,'forbackcorr',forbackcorr,'forbackWOlatencycorr',forbackWOlatencycorr,'forinvertedbackWOlatencycorr',forinvertedbackWOlatencycorr,'forwardinvertedcorr',forwardinvertedbackcorr,'latencyFor',latencyFor,'latencyBack',latencyBack);
        
    end
    
    filepath = '/Volumes/Data/DATA/Forward-Backward/Jerome/cross-correlogram';
    
    if qual_cross == 1
        
        method = 'correlation';
        
        filename = 'xcorr';
        
        save(strcat(filepath,'/xcorr/',name,'-cross-',method),'correlation');  
    
    elseif qual_cross == 2
        
        method = 'correlogram';
        
        filename = 'correlogram';
        
        save(strcat(filepath,'/correlogram/',name,'-cross-',method),'correlation');
        
    elseif qual_cross == 3
        
        method = 'correlation';
        
        filename = 'getCorrelation';
        
        save(strcat(filepath,'/getCorrelation/',name,'-cross-',method),'correlation'); 
        
    end
    
    if qual_cross == 1
        
        axes = max([length(autoforcorr) length(autobackcorr) length(forbackcorr) length(forbackWOlatencycorr) length(forinvertedbackWOlatencycorr) length(forwardinvertedbackcorr)]);
    
    elseif qual_cross == 2
        
        axes = max([length(forbackcorr) length(forbackWOlatencycorr) length(forinvertedbackWOlatencycorr) length(forwardinvertedbackcorr)]);
        
    elseif qual_cross == 3
        
        axes = max([length(autoforcorr) length(autobackcorr) length(forbackcorr) length(forbackWOlatencycorr) length(forinvertedbackWOlatencycorr) length(forwardinvertedbackcorr)]);
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = figure;
    
    if (qual_cross == 1) || (qual_cross == 3)
        
        plot(((-length(autoforcorr)/2) + 1):(length(autoforcorr)/2),autoforcorr,'k');
        hold on;
        
        plot(((-length(autobackcorr)/2) + 1):(length(autobackcorr)/2),autobackcorr,'m');
        hold on;
        
    end
    
    plot(((-length(forbackcorr)/2) + 1):(length(forbackcorr)/2),forbackcorr,'b');
    hold on;
    plot(((-length(forbackWOlatencycorr)/2) + 1):(length(forbackWOlatencycorr)/2),forbackWOlatencycorr,'r');
    hold on;
    plot(((-length(forinvertedbackWOlatencycorr)/2) + 1):(length(forinvertedbackWOlatencycorr)/2),forinvertedbackWOlatencycorr,'g');
    hold on;
    plot(((-length(forwardinvertedbackcorr)/2) + 1):(length(forwardinvertedbackcorr)/2),forwardinvertedbackcorr,'y');   
    title(strcat(name,'-latencyFor-',num2str(latencyFor),'-latencyBack-',num2str(latencyBack)));
    
    if (qual_cross == 1) || (qual_cross == 3)
        
        legend('AutoForCorr','AutoBackCorr','ForBackCorr','ForBackCorrWOLatency','ForInvertedBackCorrWOLatency','ForInvertedBackCorr');
       
    elseif qual_cross == 2
        
        legend('ForBackCorr','ForBackCorrWOLatency','ForInvertedBackCorrWOLatency','ForInvertedBackCorr');
   
    end
    
    xlim([-axes axes]);
    
    if qual_cross == 1
        
        print(f,'-depsc',strcat(filepath,'/xcorr/',name,'-',filename,'-latencyFor-',num2str(latencyFor),'-latencyBack-',num2str(latencyBack),'.eps'));
        
    elseif qual_cross == 2
    
        print(f,'-depsc',strcat(filepath,'/correlogram/',name,'-',filename,'-latencyFor-',num2str(latencyFor),'-latencyBack-',num2str(latencyBack),'.eps'));
        
    elseif qual_cross == 3
        
        print(f,'-depsc',strcat(filepath,'/getCorrelation/',name,'-',filename,'-latencyFor-',num2str(latencyFor),'-latencyBack-',num2str(latencyBack),'.eps'));
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = figure;
    
    if (qual_cross == 1) || (qual_cross == 3)
        
        plot(((-length(autoforcorr)/2) + 1):(length(autoforcorr)/2),autoforcorr,'k');
        hold on;
        
        plot(((-length(autobackcorr)/2) + 1):(length(autobackcorr)/2),autobackcorr,'m');
        hold on;
        
    end
  
    plot(((-length(forbackWOlatencycorr)/2) + 1):(length(forbackWOlatencycorr)/2),forbackWOlatencycorr,'r');
    hold on;
    plot(((-length(forinvertedbackWOlatencycorr)/2) + 1):(length(forinvertedbackWOlatencycorr)/2),forinvertedbackWOlatencycorr,'g');
   
    title(strcat(name,'-latencyFor-',num2str(latencyFor),'-latencyBack-',num2str(latencyBack)));
    
    if (qual_cross == 1) || (qual_cross == 3)
        
        legend(strcat('AutoForCorr:',num2str(max(autoforcorr)),'/','Lag:',num2str(find(autoforcorr==max(autoforcorr))),strcat('AutoBackCorr:',num2str(max(autobackcorr)),'/','Lag:',num2str(find(autobackcorr==max(autobackcorr))),strcat('ForBackCorrWOLatency',num2str(max(forbackWOlatencycorr)),'/','Lag:',num2str(find(forbackWOlatencycorr==max(forbackWOlatencycorr)))))),strcat('ForInvertedBackCorrWOLatency',num2str(max(forinvertedbackWOlatencycorr)),'/','Lag:',num2str(find(forinvertedbackWOlatencycorr==max(forinvertedbackWOlatencycorr)))));
        
    elseif qual_cross == 2
        
        legend(strcat('ForBackCorrWOLatency',num2str(max(forbackWOlatencycorr)),'/','Lag:',num2str(find(forbackWOlatencycorr==max(forbackWOlatencycorr))),strcat('ForInvertedBackCorrWOLatency',num2str(max(forinvertedbackWOlatencycorr))),'/','Lag:',num2str(find(forinvertedbackWOlatencycorr==max(forinvertedbackWOlatencycorr)))));
        
    end
    
    xlim([-axes axes]);
    
    if qual_cross == 1
        
        print(g,'-depsc',strcat(filepath,'/xcorr/',name,'-',filename,'-ForBackWOlatency-Before-After-Inverter.eps'));
        
    elseif qual_cross == 2
    
        print(g,'-depsc',strcat(filepath,'/correlogram/',name,'-',filename,'-ForBackWOlatency-Before-After-Inverter.eps'));
        
    elseif qual_cross == 3
        
        print(g,'-depsc',strcat(filepath,'/getCorrelation/',name,'-',filename,'-ForBackWOlatency-Before-After-Inverter.eps'));
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    clear correlation;
    
    clear forwardtrials;
    clear backwardtrials;
    clear autoforcorr;
    clear forbackcorr;
    clear forbackWOlatencycorr;
    clear forinvertedbackWOlatencycorr;
    clear inverted;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


