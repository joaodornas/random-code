function claudsvictor2001

disp('BEGIN');

tic

%'ice043b03_1b'
registros = { 'ice071a03_1b' };

category = { '45' '90' '135' '180' '225' '270' '315' '360' };

start_time = 0/1000;

prebaseline = 500/1000;

end_time = prebaseline + 1000/1000;

posbaseline = end_time + 2000/1000;

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/_outros/Clauds/victor2001/');

for c=1:length(registros)

    %formatSpikeData(registros,c,start_time,prebaseline,end_time,posbaseline,filepath);
    
    formatMetricSpaceGraph3(registros,c,start_time,prebaseline,end_time,posbaseline,filepath,category);
    
    formatMetricSpaceGraph5(registros,c,start_time,prebaseline,end_time,posbaseline,filepath,category);
    
    runMetricSpaceGraph3(registros,c,start_time,prebaseline,end_time,posbaseline,filepath);
    
    runMetricSpaceGraph5(registros,c,start_time,prebaseline,end_time,posbaseline,filepath);
    
end


    function formatSpikeData(registros,c,start_time,prebaseline,end_time,posbaseline,filepath)


           %%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(strcat('Load Spass...cell:',char(registros{c})));

                Spass = load(strcat('psth_',char(registros{c}),'_convolved','.mat'));


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Load Conditions Trials Labels...');

                nConditions = length(Spass.stimIds);

                for i=1:nConditions

                    trials_label(i).label = find(Spass.stimIds == i); 

                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

            %%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Set Resolution...');

                spike_times = Spass.spike_times;

                spike_times = spike_times./ 32000;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Begin Conditions Loop...');

                o = 0;
                a = 1;
                for i=1:nConditions

                    nSpikeTrains = length(spike_times(trials_label(i).label));

                    o = o + 1;

                    for j=1:nSpikeTrains

                        spiketrain = spike_times(trials_label(i).label(j),:);

                        spiketrain = spiketrain(spiketrain>0);

                        spiketrain = spiketrain(spiketrain>start_time & spiketrain <posbaseline);

                        contrast(a).orientation(o).fulltrial(j).train = spiketrain;

                        clear spiketrain;

                    end

                    if mod(o,8) == 0

                        o = 0;
                        a = a + 1; 

                    end

                end

            %%% SEPARATE SPIKETRAINS IN PARTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            disp('Separate SpikeTrains in Parts...');

            nA = length(contrast);
            for A=1:nA

                nO = length(contrast(A).orientation);

                for O=1:nO

                    nSpikeTrains = length(contrast(A).orientation(O).fulltrial);

                    spikes_vector = [];
                    for j=1:nSpikeTrains

                       spikes_vector = [spikes_vector, contrast(A).orientation(O).fulltrial(j).train];

                    end

                    spikes_vector = sort(spikes_vector);

                    contrast(A).orientation(O).spikes_vector = spikes_vector;

                    time_points = linspace(start_time, posbaseline, ( posbaseline - start_time ) * 1000);

                    [kernel.density, kernel.timePoints, kernel.optimalBinWidth, kernel.WBinsTested, kernel.Cost, kernel.confb95] = ssvkernel(spikes_vector,time_points);

                    contrast(A).orientation(O).kernel = kernel;

                    %%%%% CALCULATE LATENCY - ON  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    timePoints = length(kernel.timePoints);
                    timeHistoryBegin = prebaseline + 20/1000 ;
                    timeHistoryEnd = prebaseline + 150/1000 + 20/1000;
                    totalTime = posbaseline - start_time;
                    timeRatio = timePoints / totalTime ;

                    Ypico = max(kernel.density(timeHistoryBegin*timeRatio:timeHistoryEnd*timeRatio));

                    for k=timeHistoryBegin*timeRatio:timeHistoryEnd*timeRatio

                        if kernel.density(k) == Ypico

                            Xpico = k;

                        end

                    end

                    ae = median(kernel.density(1:(timeHistoryBegin - 20/1000)*timeRatio));

                    YM = ( (Ypico - ae) / 2 ) + ae ;
                    kernel.density = kernel.density';
                    V = kernel.density(timeHistoryBegin*timeRatio:Xpico);
                    XM = dsearchn(V,YM);
                    %XM = XM + timeHistoryBegin*timeRatio;
                    latencyon = kernel.timePoints(XM);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    contrast(A).orientation(O).latencyon = latencyon;

                    %%%%% CALCULATE LATENCY  - OFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Ypico = max(kernel.density((prebaseline+latencyon)*timeRatio:(prebaseline+500/1000)*timeRatio));

                    latencyON = round(latencyon*1000);

                    for k=((prebaseline*timeRatio)+latencyON):(prebaseline+500/1000)*timeRatio

                        if kernel.density(k) == Ypico

                            Xpico = k;

                        end

                    end

                    bl = median(kernel.density((prebaseline+500/1000)*timeRatio:end_time*timeRatio));

                    YM = ( (Ypico - bl) / 2 ) + bl ;
                    %kernel.density = kernel.density';
                    V = kernel.density(Xpico:(prebaseline+500/1000)*timeRatio);
                    XM = dsearchn(V,YM);
                    %XM = XM + timeHistoryBegin*timeRatio;
                    latencyoff = kernel.timePoints(XM);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    contrast(A).orientation(O).latencyoff = (latencyoff + Xpico) / 1000;

                    for j=1:nSpikeTrains

                       spiketrain = contrast(A).orientation(O).fulltrial(j).train;

                       fullresponse = spiketrain(spiketrain>=prebaseline & spiketrain<=end_time);

                       contrast(A).orientation(O).fullresponse(j).train = fullresponse;

                       wolatency = fullresponse - contrast(A).orientation(O).latencyon;

                       wolatency = wolatency(wolatency>=prebaseline & wolatency<=end_time);

                       contrast(A).orientation(O).wolatency(j).train = wolatency;

                       if length(wolatency) < 1

                           contrast(A).orientation(O).spikealone(j).train = [];   

                       else

                           contrast(A).orientation(O).spikealone(j).train = wolatency(1); 

                       end
                       
                       transient = spiketrain(spiketrain>=contrast(A).orientation(O).latencyon & spiketrain<=contrast(A).orientation(O).latencyoff);
                       
                       transient = transient - contrast(A).orientation(O).latencyon;
                       
                       transient = transient(transient>=prebaseline & transient<=end_time);
                       
                       contrast(A).orientation(O).transient(j).train = transient;
                       
                       tonic = spiketrain(spiketrain>=contrast(A).orientation(O).latencyoff & spiketrain<=end_time);
                       
                       tonic = tonic - contrast(A).orientation(O).latencyon;
                       
                       tonic = tonic(tonic>=contrast(A).orientation(O).latencyoff & tonic<=end_time);
                       
                       contrast(A).orientation(O).tonic(j).train = tonic;
                       
                       offtrain = spiketrain(spiketrain>=end_time & spiketrain<=posbaseline);
                       
                       contrast(A).orientation(O).off(j).train = offtrain;

                       clear spiketrain;
                       clear fullresponse;
                       clear transient;
                       clear offtrain;
                       clear tonic;

                    end

                end

            end

            save(strcat(filepath,char(registros{c}),'-contrast'),'contrast');   
    end 
    

    function formatMetricSpaceGraph3(registros,c,start_time,prebaseline,end_time,posbaseline,filepath,category)
    
        pack = load(strcat(char(registros{c}),'-contrast.mat'));

        nA = length(pack.contrast);
        for A=1:nA

            %%%   SET METRIC SPACE PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Set Metric Space parameters...');

            %Escala de tempo em segundos
            time_scale = 1;

            %Escala da resolu??o temporal em segundos
            time_resolution = 0.0001;

            %Sistema Internacional de Medidas
            si_prefix = 1;

            %N?mero de classes de est?mulos
            M = length(pack.contrast(A).orientation);

            %N?mero de s?tios (canais)
            N = 1;

            %N?mero de trials (repeti??es por condi??o)
            P = 10;

            %Tempos inicial e final do trial
            %start_time = 0;
            %end_time = 0;

            %N?mero de pontos no vetor list
            %Q = 0;

            categoriesFull(1:M) = struct('label','','P',P,'trials',zeros(P));
            categoriesWOLatency(1:M) = struct('label','','P',P,'trials',zeros(P));
            categoriesSpikeAlone(1:M) = struct('label','','P',P,'trials',zeros(P));

            nO = length(pack.contrast(A).orientation);

            for O=1:nO

                nSpikeTrain = length(pack.contrast(A).orientation(O).fulltrial);

                %trialsFull(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);
                %trialsWOLatency(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);
                %trialsSpikeAlone(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);

                F = 1;
                W = 1;
                S = 1;
                for j=1:nSpikeTrain

                    spikesfull = pack.contrast(A).orientation(O).fullresponse(j).train;
                    spikeswolatency = pack.contrast(A).orientation(O).wolatency(j).train;
                    spikesspikealone = pack.contrast(A).orientation(O).spikealone(j).train;


                    if size(spikesfull,2)>0

                        trialsFull(F).start_time = prebaseline;
                        trialsFull(F).end_time = end_time;
                        trialsFull(F).Q = size(spikesfull,2);
                        trialsFull(F).list = spikesfull;
                        
                        F = F + 1;

                    end

                    if size(spikeswolatency,2)>0

                        trialsWOLatency(W).start_time = prebaseline;
                        trialsWOLatency(W).end_time = end_time;
                        trialsWOLatency(W).Q = size(spikeswolatency,2);
                        trialsWOLatency(W).list = spikeswolatency;

                        W = W + 1;

                    end

                    if size(spikesspikealone,2)>0

                        trialsSpikeAlone(S).start_time = prebaseline;
                        trialsSpikeAlone(S).end_time = end_time;
                        trialsSpikeAlone(S).Q = size(spikesspikealone,2);
                        trialsSpikeAlone(S).list = spikesspikealone;    

                        S = S + 1;

                    end

                end

            trialsFull = trialsFull';
            trialsWOLatency = trialsWOLatency';
            trialsSpikeAlone = trialsSpikeAlone';

            categoriesFull(O).label = category{O};
            categoriesFull(O).P = F-1;
            categoriesFull(O).trials = trialsFull;
            
            categoriesWOLatency(O).label = category{O};
            categoriesWOLatency(O).P = W-1;
            categoriesWOLatency(O).trials = trialsWOLatency;

            categoriesSpikeAlone(O).label = category{O};
            categoriesSpikeAlone(O).P = S-1;
            categoriesSpikeAlone(O).trials = trialsSpikeAlone;

            categoriesFull = categoriesFull';
            categoriesWOLatency = categoriesWOLatency';
            categoriesSpikeAlone = categoriesSpikeAlone';

            clear trialsFull;
            clear trialsWOLatency;
            clear trialsSpikeAlone;

            end

            sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

            XFull = struct('M',M,'N',N,'sites',sites,'categories',categoriesFull);
            XWOLatency = struct('M',M,'N',N,'sites',sites,'categories',categoriesWOLatency);
            XSpikeAlone = struct('M',M,'N',N,'sites',sites,'categories',categoriesSpikeAlone);

            fileFull{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Full-Response');
            fileWOLatency{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','WOLatency');
            fileSpikeAlone{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','SpikeAlone');

            stawrite(XFull,char(fileFull{A}));
            stawrite(XWOLatency,char(fileWOLatency{A}));
            stawrite(XSpikeAlone,char(fileSpikeAlone{A}));

            clear categoriesFull;
            clear categoriesWOLatency;
            clear categoriesSpikeAlone;

        end
        
    end


    function formatMetricSpaceGraph5(registros,c,start_time,prebaseline,end_time,posbaseline,filepath,category)

            pack = load(strcat(char(registros{c}),'-contrast.mat'));

            nA = length(pack.contrast);
            for A=1:nA

                %%%   SET METRIC SPACE PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Set Metric Space parameters...');

                %Escala de tempo em segundos
                time_scale = 1;

                %Escala da resolu??o temporal em segundos
                time_resolution = 0.0001;

                %Sistema Internacional de Medidas
                si_prefix = 1;

                %N?mero de classes de est?mulos
                M = length(pack.contrast(A).orientation);

                %N?mero de s?tios (canais)
                N = 1;

                %N?mero de trials (repeti??es por condi??o)
                P = 10;

                %Tempos inicial e final do trial
                %start_time = 0;
                %end_time = 0;

                %N?mero de pontos no vetor list
                %Q = 0;

                %categoriesTransient(1:M) = struct('label','','P',P,'trials',zeros(P));
                %categoriesTonic(1:M) = struct('label','','P',P,'trials',zeros(P));
                %categoriesOff(1:M) = struct('label','','P',P,'trials',zeros(P));

                nO = length(pack.contrast(A).orientation);

                R = 1;
                for O=1:nO

                    nSpikeTrain = length(pack.contrast(A).orientation(O).fulltrial);

                    %trialsTransient(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);
                    %trialsTonic(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);
                    %trialsOff(1:nSpikeTrain) = struct('start_time',prebaseline,'end_time',end_time,'Q',0,'list',[]);

                    F = 1;
                    W = 1;
                    S = 1;
                    for j=1:nSpikeTrain

                        spikestransient = pack.contrast(A).orientation(O).transient(j).train;
                        spikestonic = pack.contrast(A).orientation(O).tonic(j).train;
                        spikesoff = pack.contrast(A).orientation(O).off(j).train;


                        if size(spikestransient,2)>0

                            trialsTransient(F).start_time = prebaseline;
                            trialsTransient(F).end_time = end_time;
                            trialsTransient(F).Q = size(spikestransient,2);
                            trialsTransient(F).list = spikestransient;

                            F = F + 1;

                        end

                        if size(spikestonic,2)>0

                            trialsTonic(W).start_time = prebaseline;
                            trialsTonic(W).end_time = end_time;
                            trialsTonic(W).Q = size(spikestonic,2);
                            trialsTonic(W).list = spikestonic;

                            W = W + 1;

                        end

                        if size(spikesoff,2)>0

                            trialsOff(S).start_time = end_time;
                            trialsOff(S).end_time = posbaseline;
                            trialsOff(S).Q = size(spikesoff,2);
                            trialsOff(S).list = spikesoff;    

                            S = S + 1;

                        end

                    end
                    
                    if (F > 1) && (W > 1) && (S > 1)
                        
                        trialsTransient = trialsTransient';
                        trialsTonic = trialsTonic';
                        trialsOff = trialsOff';

                        categoriesTransient(R).label = category{R};
                        categoriesTransient(R).P = F-1;
                        categoriesTransient(R).trials = trialsTransient;

                        categoriesTonic(R).label = category{R};
                        categoriesTonic(R).P = W-1;
                        categoriesTonic(R).trials = trialsTonic;

                        categoriesOff(R).label = category{R};
                        categoriesOff(R).P = S-1;
                        categoriesOff(R).trials = trialsOff;

                        categoriesTransient = categoriesTransient';
                        categoriesTonic = categoriesTonic';
                        categoriesOff = categoriesOff';

                        clear trialsTransient;
                        clear trialsTonic;
                        clear trialsOff;
                        
                        R = R + 1;
                        
                    end
                  
                end

                sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

                XTransient = struct('M',R-1,'N',N,'sites',sites,'categories',categoriesTransient);
                XTonic = struct('M',R-1,'N',N,'sites',sites,'categories',categoriesTonic);
                XOff = struct('M',R-1,'N',N,'sites',sites,'categories',categoriesOff);

                fileTransient{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Transient');
                fileTonic{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Tonic');
                fileOff{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Off');
                
                stawrite(XTransient,char(fileTransient{A}));
                stawrite(XTonic,char(fileTonic{A}));
                stawrite(XOff,char(fileOff{A}));

                clear categoriesTransient;
                clear categoriesTonic;
                clear categoriesOff;

            end

      end


    function runMetricSpaceGraph3(registros,c,start_time,prebaseline,end_time,posbaseline,filepath)

        pack = load(strcat(char(registros{c}),'-contrast.mat'));

        nA = length(pack.contrast);
        for A=1:nA

            P = 10;
            
            disp('READING DATA');
            
            fileFull{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Full-Response');
            fileWOLatency{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','WOLatency');
            fileSpikeAlone{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','SpikeAlone');
            
            % READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
            XFull = staread(strrep(strcat(char(fileFull{A}),'.stam'),'/',filesep));
            XWOLatency = staread(strrep(strcat(char(fileWOLatency{A}),'.stam'),'/',filesep));
            XSpikeAlone = staread(strrep(strcat(char(fileSpikeAlone{A}),'.stam'),'/',filesep));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp('SETTING OPTIONS');
            % OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            opts.entropy_estimation_method = {'plugin','tpmc','jack'};
            %opts.variance_estimation_method = {'jack'};

            %opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
            opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
            %opts.unoccupied_bins_strategy = 1; % Use all bins

            opts.parallel = 1;
            opts.possible_words = 'unique';

            opts.start_time = start_time;
            opts.end_time = end_time;
            opts.shift_cost = [0 2.^(-4:9)];
            %opts.label_cost = [0 1 2];
            opts.clustering_exponent = -2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            S = P*10;

            disp('METRIC TIMING');
            % METRIC TIMING + SHUFFLE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            opts.metric_family = 0;

            disp('...computing metric space timing');
            [out_TFull,opts_used] = metric(XFull,opts);
            [out_TWOLatency,opts_used] = metric(XWOLatency,opts);
            [out_TSpikeAlone,opts_used] = metric(XSpikeAlone,opts);

            for q_idx=1:length(opts.shift_cost)

                info_plugin_TFull(q_idx) = out_TFull(q_idx).table.information(1).value;
                info_tpmc_TFull(q_idx) = out_TFull(q_idx).table.information(2).value;
                info_jack_TFull(q_idx) = out_TFull(q_idx).table.information(3).value;

                info_plugin_TWOLatency(q_idx) = out_TWOLatency(q_idx).table.information(1).value;
                info_tpmc_TWOLatency(q_idx) = out_TWOLatency(q_idx).table.information(2).value;
                info_jack_TWOLatency(q_idx) = out_TWOLatency(q_idx).table.information(3).value;

                info_plugin_TSpikeAlone(q_idx) = out_TSpikeAlone(q_idx).table.information(1).value;
                info_tpmc_TSpikeAlone(q_idx) = out_TSpikeAlone(q_idx).table.information(2).value;
                info_jack_TSpikeAlone(q_idx) = out_TSpikeAlone(q_idx).table.information(3).value;

            end

            [max_info_plugin_TFull,max_info_plugin_idx_TFull] = max(info_plugin_TFull);
            [max_info_plugin_TWOLatency,max_info_plugin_idx_TWOLatency] = max(info_plugin_TWOLatency);
            [max_info_plugin_TSpikeAlone,max_info_plugin_idx_TSpikeAlone] = max(info_plugin_TSpikeAlone);

            [max_info_tpmc_TFull,max_info_tpmc_idx_TFull] = max(info_tpmc_TFull);
            [max_info_tpmc_TWOLatency,max_info_tpmc_idx_TWOLatency] = max(info_tpmc_TWOLatency);
            [max_info_tpmc_TSpikeAlone,max_info_tpmc_idx_TSpikeAlone] = max(info_tpmc_TSpikeAlone);

            [max_info_jack_TFull,max_info_jack_idx_TFull] = max(info_jack_TFull);
            [max_info_jack_TWOLatency,max_info_jack_idx_TWOLatency] = max(info_jack_TWOLatency);
            [max_info_jack_TSpikeAlone,max_info_jack_idx_TSpikeAlone] = max(info_jack_TSpikeAlone);

            Hcount_plugin_TFull = info_plugin_TFull(1);
            Hcount_tpmc_TFull = info_tpmc_TFull(1);
            Hcount_jack_TFull = info_jack_TFull(1);

            Hcount_plugin_TWOLatency = info_plugin_TWOLatency(1);
            Hcount_tpmc_TWOLatency = info_tpmc_TWOLatency(1);
            Hcount_jack_TWOLatency = info_jack_TWOLatency(1);

            Hcount_plugin_TSpikeAlone = info_plugin_TSpikeAlone(1);
            Hcount_tpmc_TSpikeAlone = info_tpmc_TSpikeAlone(1);
            Hcount_jack_TSpikeAlone = info_jack_TSpikeAlone(1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp('CALCULATE H-BIAS VIA SHUFFLE');
            % CALCULATE H-BIAS VIA SHUFFLE %%%%%%%%%%%%%%%

            disp('...computing metric space timing shuffle');

            opts.entropy_estimation_method = {'plugin'};

            rand('state',0);
            opts.metric_family = 0;
            [Y_TFull,SHUF_TFull,opts_used] = metric_shuf(XFull,opts,S);
            [Y_TWOLatency,SHUF_TWOLatency,opts_used] = metric_shuf(XWOLatency,opts,S);
            [Y_TSpikeAlone,SHUF_TSpikeAlone,opts_used] = metric_shuf(XSpikeAlone,opts,S);

            SHUF_TFull = SHUF_TFull';
            SHUF_TWOLatency = SHUF_TWOLatency';
            SHUF_TSpikeAlone = SHUF_TSpikeAlone';

            for q_idx=1:length(opts.shift_cost)

                for i=1:S

                    HShuffle_TFull(i,q_idx) = SHUF_TFull(i,q_idx).table.information.value;
                    HShuffle_TWOLatency(i,q_idx) = SHUF_TWOLatency(i,q_idx).table.information.value;
                    HShuffle_TSpikeAlone(i,q_idx) = SHUF_TSpikeAlone(i,q_idx).table.information.value;

                end

            end

            HBias_TFull = mean(HShuffle_TFull,1);
            HBias_std_TFull = std(HShuffle_TFull,[],1);

            HBias_TWOLatency = mean(HShuffle_TWOLatency,1);
            HBias_std_TWOLatency = std(HShuffle_TWOLatency,[],1);

            HBias_TSpikeAlone = mean(HShuffle_TSpikeAlone,1);
            HBias_std_TSpikeAlone = std(HShuffle_TSpikeAlone,[],1);

            %%% leave-one-out Jackknife 
            disp('leave-one-out Jackknife');

            opts.entropy_estimation_method = {'plugin'};
            opts.metric_family = 0;

            [out_unjkFull,jkFull,opts_used] = metric_jack(XFull,opts);
            [out_unjkWOLatency,jkWOLatency,opts_used] = metric_jack(XWOLatency,opts);
            [out_unjkSpikeAlone,jkSpikeAlone,opts_used] = metric_jack(XSpikeAlone,opts);

            P_totalFull = size(jkFull,1);
            P_totalWOLatency = size(jkWOLatency,1);
            P_totalSpikeAlone = size(jkSpikeAlone,1);

            temp_info_jkFull = zeros(P_totalFull,length(opts.shift_cost));
            temp_info_jkWOLatency = zeros(P_totalWOLatency,length(opts.shift_cost));
            temp_info_jkSpikeAlone = zeros(P_totalSpikeAlone,length(opts.shift_cost));

            for q_idx=1:length(opts.shift_cost)

              info_unjkFull(q_idx)= out_unjkFull(q_idx).table.information.value;
              info_unjkWOLatency(q_idx)= out_unjkWOLatency(q_idx).table.information.value;
              info_unjkSpikeAlone(q_idx)= out_unjkSpikeAlone(q_idx).table.information.value;

              for p=1:P_totalFull

                temp_info_jkFull(p,q_idx) = jkFull(p,q_idx).table.information.value;

              end

              for p=1:P_totalWOLatency

                temp_info_jkWOLatency(p,q_idx) = jkWOLatency(p,q_idx).table.information.value;

              end

              for p=1:P_totalSpikeAlone

                temp_info_jkSpikeAlone(p,q_idx) = jkSpikeAlone(p,q_idx).table.information.value;

              end

            end

            info_jk_semFull = sqrt((P_totalFull-1)*var(temp_info_jkFull,1,1));
            info_jk_semWOLatency = sqrt((P_totalWOLatency-1)*var(temp_info_jkWOLatency,1,1));
            info_jk_semSpikeAlone = sqrt((P_totalSpikeAlone-1)*var(temp_info_jkSpikeAlone,1,1)); 
            
            x = [3;3;3;3;3];
            
            q = opts.shift_cost.';
            
            options.MaxFunEvals = 500000;
            options.MaxIter = 500000;
            options.TolFun = 1.000000e-06;
            options.TolX = 1.000000e-06;
            
            [HfitFull resnormFull] = lsqcurvefit(@Hfit,x,q,info_tpmc_TFull.',[],[],options);
            [HfitWOLatency resnormWOLatency] = lsqcurvefit(@Hfit,x,q,info_tpmc_TWOLatency.',[],[],options);
            [HfitSpikeAlone resnormSpikeAlone] = lsqcurvefit(@Hfit,x,q,info_tpmc_TSpikeAlone.',[],[],options);
            
            HfitFullhalf = max_info_tpmc_TFull / 2;
            HfitWOLatencyhalf = max_info_tpmc_TWOLatency / 2;
            HfitSpikeAlonehalf = max_info_tpmc_TSpikeAlone / 2;            
            
            for y=1:512
                
                HF(y) = HfitFull(1)*((1 + HfitFull(2)*y.^HfitFull(3))./(1 + HfitFull(4)*y.^HfitFull(5)));
                HWOL(y) = HfitWOLatency(1)*((1 + HfitWOLatency(2)*y.^HfitWOLatency(3))./(1 + HfitWOLatency(4)*y.^HfitWOLatency(5)));
                HSA(y) = HfitSpikeAlone(1)*((1 + HfitSpikeAlone(2)*y.^HfitSpikeAlone(3))./(1 + HfitSpikeAlone(4)*y.^HfitSpikeAlone(5)));  
    
            end
            
            [d qFull] = min(abs(HF - HfitFullhalf));
            [d qWOLatency] = min(abs(HWOL - HfitWOLatencyhalf));
            [d qSpikeAlone] = min(abs(HSA - HfitSpikeAlonehalf));
                         
            tempPrecisionFull = 2000/qFull;
            tempPrecisionWOLatency = 2000/qWOLatency;
            tempPrecisionSpikeAlone = 2000/qSpikeAlone;
            
            percTempFull = 100 * (max_info_tpmc_TFull - Hcount_plugin_TFull) / max_info_tpmc_TFull;
            percTempWOLatency = 100 * (max_info_tpmc_TWOLatency - Hcount_tpmc_TWOLatency) / max_info_tpmc_TWOLatency;
            percTempSpikeAlone = 100 * (max_info_tpmc_TSpikeAlone - Hcount_tpmc_TSpikeAlone) / max_info_tpmc_TSpikeAlone;
            
            max_info = struct('max_info_plugin_idx_TFull',max_info_plugin_idx_TFull,'max_info_plugin_TFull',max_info_plugin_TFull,'max_info_plugin_idx_TWOLatency',max_info_plugin_idx_TWOLatency,'max_info_plugin_TWOLatency',max_info_plugin_TWOLatency,'max_info_plugin_idx_TSpikeAlone',max_info_plugin_idx_TSpikeAlone,'max_info_plugin_TSpikeAlone',max_info_plugin_TSpikeAlone,'Hcount_plugin_info_TFull',Hcount_plugin_TFull,'Hcount_plugin_info_TWOLatency',Hcount_plugin_TWOLatency,'Hcount_plugin_info_TSpikeAlone',Hcount_plugin_TSpikeAlone,'max_info_tpmc_idx_TFull',max_info_tpmc_idx_TFull,'max_info_tpmc_TFull',max_info_tpmc_TFull,'max_info_tpmc_idx_TWOLatency',max_info_tpmc_idx_TWOLatency,'max_info_tpmc_TWOLatency',max_info_tpmc_TWOLatency,'max_info_tpmc_idx_TSpikeAlone',max_info_tpmc_idx_TSpikeAlone,'max_info_tpmc_TSpikeAlone',max_info_tpmc_TSpikeAlone,'Hcount_tpmc_info_TFull',Hcount_plugin_TFull,'Hcount_tpmc_info_TWOLatency',Hcount_tpmc_TWOLatency,'Hcount_tpmc_info_TSpikeAlone',Hcount_tpmc_TSpikeAlone,'max_info_jack_idx_TFull',max_info_jack_idx_TFull,'max_info_jack_TFull',max_info_jack_TFull,'Hcount_jack_info_TFull',Hcount_jack_TFull,'HBias_TFull',HBias_TFull,'HBias_std_TFull',HBias_std_TFull,'max_info_jack_idx_TWOLatency',max_info_jack_idx_TWOLatency,'max_info_jack_TTonic',max_info_jack_TWOLatency,'Hcount_jack_info_TWOLatency',Hcount_jack_TWOLatency,'HBias_TWOLatency',HBias_TWOLatency,'HBias_std_TWOLatency',HBias_std_TWOLatency,'max_info_jack_idx_TSpikeAlone',max_info_jack_idx_TSpikeAlone,'max_info_jack_TSpikeAlone',max_info_jack_TSpikeAlone,'Hcount_jack_info_TSpikeAlone',Hcount_jack_TSpikeAlone,'HBias_TSpikeAlone',HBias_TSpikeAlone,'HBias_std_TSpikeAlone',HBias_std_TSpikeAlone,'info_unjkFull',info_unjkFull,'info_unjkWOLatency',info_unjkWOLatency,'info_unjkSpikeAlone',info_unjkSpikeAlone,'HfitFull',HfitFull,'HfitWOLatency',HfitWOLatency,'HfitSpikeAlone',HfitSpikeAlone,'qFull',qFull,'qWOLatency',qWOLatency,'qSpikeAlone',qSpikeAlone,'tempPrecisionFull',tempPrecisionFull,'tempPrecisionWOLatency',tempPrecisionWOLatency,'tempPrecisionSpikeAlone',tempPrecisionSpikeAlone,'percTempFull',percTempFull,'percTempWOLatency',percTempWOLatency,'percTempSpikeAlone',percTempSpikeAlone);

            metric_analysis = struct('XFull', XFull, 'XWOLatency', XWOLatency, 'XSpikeAlone', XSpikeAlone, 'out_TFull', out_TFull, 'out_TWOLatency', out_TWOLatency, 'out_TSpikeAlone', out_TSpikeAlone, 'max_info', max_info, 'opts', opts, 'Y_TFull', Y_TFull, 'Y_TWOLatency', Y_TWOLatency, 'Y_TSpikeAlone', Y_TSpikeAlone, 'SHUF_TFull', SHUF_TFull, 'SHUF_TWOLatency', SHUF_TWOLatency,'SHUF_TSpikeAlone', SHUF_TSpikeAlone);

            save(strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-metric_analysis-graph-3'),'metric_analysis');  

            disp('PLOTTING INFO TIMING-INTERVAL');
            % INFO TIMING-INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f(A) = figure;
            %plot(1:length(opts.shift_cost),info_tpmc_TTransient,'bo');
            errorbar(1:length(opts.shift_cost),info_tpmc_TFull,2*info_jk_semFull,'bo');
            hold on;
            
            %plot(1:length(opts.shift_cost),info_tpmc_TTonic,'ro');
            errorbar(1:length(opts.shift_cost),info_tpmc_TWOLatency,2*info_jk_semWOLatency,'ro');
            
            %plot(1:length(opts.shift_cost),info_tpmc_TOff,'go');
            errorbar(1:length(opts.shift_cost),info_tpmc_TSpikeAlone,2*info_jk_semSpikeAlone,'go');
            hold off;

            set(gca,'xtick',1:length(opts.shift_cost));
            set(gca,'xticklabel',opts.shift_cost);
            set(gca,'xlim',[1 length(opts.shift_cost)]);
            set(gca,'ylim',[-0.5 2.5]);

            xlabel('Temporal precision (1/sec)');
            ylabel('Information (bits)');

            legend('Full Response','WO Latency','1st Spike Alone');

            print(f(A),'-depsc',strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-plot-info-timing-latency'));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end     
        
    end


    function runMetricSpaceGraph5(registros,c,start_time,prebaseline,end_time,posbaseline,filepath)

        pack = load(strcat(char(registros{c}),'-contrast.mat'));

        nA = length(pack.contrast);
        for A=1:nA

            P = 10;
            
            disp('READING DATA');
            
            fileTransient{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Transient');
            fileTonic{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Tonic');
            fileOff{A} = strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-','Off');
            
            % READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
            XTransient = staread(strrep(strcat(char(fileTransient{A}),'.stam'),'/',filesep));
            XTonic = staread(strrep(strcat(char(fileTonic{A}),'.stam'),'/',filesep));
            XOff = staread(strrep(strcat(char(fileOff{A}),'.stam'),'/',filesep));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp('SETTING OPTIONS');
            % OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            opts.entropy_estimation_method = {'plugin','tpmc','jack'};
            %opts.variance_estimation_method = {'jack'};

            %opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
            opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
            %opts.unoccupied_bins_strategy = 1; % Use all bins

            opts.parallel = 1;
            opts.possible_words = 'unique';

            opts.start_time = start_time;
            opts.end_time = end_time;
            opts.shift_cost = [0 2.^(-4:9)];
            %opts.label_cost = [0 1 2];
            opts.clustering_exponent = -2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            S = P*10;

            disp('METRIC TIMING');
            % METRIC TIMING + SHUFFLE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            opts.metric_family = 0;

            disp('...computing metric space timing');
            [out_TTransient,opts_used] = metric(XTransient,opts);
            [out_TTonic,opts_used] = metric(XTonic,opts);
            [out_TOff,opts_used] = metric(XOff,opts);

            for q_idx=1:length(opts.shift_cost)

                info_plugin_TTransient(q_idx) = out_TTransient(q_idx).table.information(1).value;
                info_tpmc_TTransient(q_idx) = out_TTransient(q_idx).table.information(2).value;
                info_jack_TTransient(q_idx) = out_TTransient(q_idx).table.information(3).value;

                info_plugin_TTonic(q_idx) = out_TTonic(q_idx).table.information(1).value;
                info_tpmc_TTonic(q_idx) = out_TTonic(q_idx).table.information(2).value;
                info_jack_TTonic(q_idx) = out_TTonic(q_idx).table.information(3).value;

                info_plugin_TOff(q_idx) = out_TOff(q_idx).table.information(1).value;
                info_tpmc_TOff(q_idx) = out_TOff(q_idx).table.information(2).value;
                info_jack_TOff(q_idx) = out_TOff(q_idx).table.information(3).value;

            end

            [max_info_plugin_TTransient,max_info_plugin_idx_TTransient] = max(info_plugin_TTransient);
            [max_info_plugin_TTonic,max_info_plugin_idx_TTonic] = max(info_plugin_TTonic);
            [max_info_plugin_TOff,max_info_plugin_idx_TOff] = max(info_plugin_TOff);

            [max_info_tpmc_TTransient,max_info_tpmc_idx_TTransient] = max(info_tpmc_TTransient);
            [max_info_tpmc_TTonic,max_info_tpmc_idx_TTonic] = max(info_tpmc_TTonic);
            [max_info_tpmc_TOff,max_info_tpmc_idx_TOff] = max(info_tpmc_TOff);

            [max_info_jack_TTransient,max_info_jack_idx_TTransient] = max(info_jack_TTransient);
            [max_info_jack_TTonic,max_info_jack_idx_TTonic] = max(info_jack_TTonic);
            [max_info_jack_TOff,max_info_jack_idx_TOff] = max(info_jack_TOff);

            Hcount_plugin_TTransient = info_plugin_TTransient(1);
            Hcount_plugin_TTransient = info_tpmc_TTransient(1);
            Hcount_jack_TTransient = info_jack_TTransient(1);

            Hcount_plugin_TTonic = info_plugin_TTonic(1);
            Hcount_tpmc_TTonic = info_tpmc_TTonic(1);
            Hcount_jack_TTonic = info_jack_TTonic(1);

            Hcount_plugin_TOff = info_plugin_TOff(1);
            Hcount_tpmc_TOff = info_tpmc_TOff(1);
            Hcount_jack_TOff = info_jack_TOff(1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp('CALCULATE H-BIAS VIA SHUFFLE');
            % CALCULATE H-BIAS VIA SHUFFLE %%%%%%%%%%%%%%%

            disp('...computing metric space timing shuffle');

            opts.entropy_estimation_method = {'plugin'};

            rand('state',0);
            opts.metric_family = 0;
            [Y_TTransient,SHUF_TTransient,opts_used] = metric_shuf(XTransient,opts,S);
            [Y_TTonic,SHUF_TTonic,opts_used] = metric_shuf(XTonic,opts,S);
            [Y_TOff,SHUF_TOff,opts_used] = metric_shuf(XOff,opts,S);

            SHUF_TTransient = SHUF_TTransient';
            SHUF_TTonic = SHUF_TTonic';
            SHUF_TOff = SHUF_TOff';

            for q_idx=1:length(opts.shift_cost)

                for i=1:S

                    HShuffle_TTransient(i,q_idx) = SHUF_TTransient(i,q_idx).table.information.value;
                    HShuffle_TTonic(i,q_idx) = SHUF_TTonic(i,q_idx).table.information.value;
                    HShuffle_TOff(i,q_idx) = SHUF_TOff(i,q_idx).table.information.value;

                end

            end

            HBias_TTransient = mean(HShuffle_TTransient,1);
            HBias_std_TTransient = std(HShuffle_TTransient,[],1);

            HBias_TTonic = mean(HShuffle_TTonic,1);
            HBias_std_TTonic = std(HShuffle_TTonic,[],1);

            HBias_TOff = mean(HShuffle_TOff,1);
            HBias_std_TOff = std(HShuffle_TOff,[],1);

            %%% leave-one-out Jackknife 
            disp('leave-one-out Jackknife');

            opts.entropy_estimation_method = {'plugin'};
            opts.metric_family = 0;

            [out_unjkTransient,jkTransient,opts_used] = metric_jack(XTransient,opts);
            [out_unjkTonic,jkTonic,opts_used] = metric_jack(XTonic,opts);
            [out_unjkOff,jkOff,opts_used] = metric_jack(XOff,opts);

            P_totalTransient = size(jkTransient,1);
            P_totalTonic = size(jkTonic,1);
            P_totalOff = size(jkOff,1);

            temp_info_jkTransient = zeros(P_totalTransient,length(opts.shift_cost));
            temp_info_jkTonic = zeros(P_totalTonic,length(opts.shift_cost));
            temp_info_jkOff = zeros(P_totalOff,length(opts.shift_cost));

            for q_idx=1:length(opts.shift_cost)

              info_unjkTransient(q_idx)= out_unjkTransient(q_idx).table.information.value;
              info_unjkTonic(q_idx)= out_unjkTonic(q_idx).table.information.value;
              info_unjkOff(q_idx)= out_unjkOff(q_idx).table.information.value;

              for p=1:P_totalTransient

                temp_info_jkTransient(p,q_idx) = jkTransient(p,q_idx).table.information.value;

              end

              for p=1:P_totalTonic

                temp_info_jkTonic(p,q_idx) = jkTonic(p,q_idx).table.information.value;

              end

              for p=1:P_totalOff

                temp_info_jkOff(p,q_idx) = jkOff(p,q_idx).table.information.value;

              end

            end

            info_jk_semTransient = sqrt((P_totalTransient-1)*var(temp_info_jkTransient,1,1));
            info_jk_semTonic = sqrt((P_totalTonic-1)*var(temp_info_jkTonic,1,1));
            info_jk_semOff = sqrt((P_totalOff-1)*var(temp_info_jkOff,1,1)); 
            
            x = [3;3;3;3;3];
            
            q = opts.shift_cost.';
            
            options.MaxFunEvals = 500000;
            options.MaxIter = 500000;
            options.TolFun = 1.000000e-06;
            options.TolX = 1.000000e-06;
            
            [HfitTransient resnormTransient] = lsqcurvefit(@Hfit,x,q,info_tpmc_TTransient.',[],[],options);
            [HfitTonic resnormTonic] = lsqcurvefit(@Hfit,x,q,info_tpmc_TTonic.',[],[],options);
            [HfitOff resnormOff] = lsqcurvefit(@Hfit,x,q,info_tpmc_TOff.',[],[],options);
            
            HfitTransienthalf = max_info_tpmc_TTransient / 2;
            HfitTonichalf = max_info_tpmc_TTonic / 2;
            HfitOffhalf = max_info_tpmc_TOff / 2;            
            
            for y=1:512
                
                HTRA(y) = HfitTransient(1)*((1 + HfitTransient(2)*y.^HfitTransient(3))./(1 + HfitTransient(4)*y.^HfitTransient(5)));
                HTO(y) = HfitTonic(1)*((1 + HfitTonic(2)*y.^HfitTonic(3))./(1 + HfitTonic(4)*y.^HfitTonic(5)));
                HOFF(y) = HfitOff(1)*((1 + HfitOff(2)*y.^HfitOff(3))./(1 + HfitOff(4)*y.^HfitOff(5)));  
    
            end
            
            [d qTransient] = min(abs(HTRA - HfitTransienthalf));
            [d qTonic] = min(abs(HTO - HfitTonichalf));
            [d qOff] = min(abs(HOFF - HfitOffhalf));
                         
            tempPrecisionTransient = 2000/qTransient;
            tempPrecisionTonic = 2000/qTonic;
            tempPrecisionOff = 2000/qOff;
            
            percTempTransient = 100 * (max_info_tpmc_TTransient - Hcount_plugin_TTransient) / max_info_tpmc_TTransient;
            percTempTonic = 100 * (max_info_tpmc_TTonic - Hcount_tpmc_TTonic) / max_info_tpmc_TTonic;
            percTempOff = 100 * (max_info_tpmc_TOff - Hcount_tpmc_TOff) / max_info_tpmc_TOff;
            
            max_info = struct('max_info_plugin_idx_TTransient',max_info_plugin_idx_TTransient,'max_info_plugin_TTransient',max_info_plugin_TTransient,'max_info_plugin_idx_TTonic',max_info_plugin_idx_TTonic,'max_info_plugin_TTonic',max_info_plugin_TTonic,'max_info_plugin_idx_TOff',max_info_plugin_idx_TOff,'max_info_plugin_TOff',max_info_plugin_TOff,'Hcount_plugin_info_TTransient',Hcount_plugin_TTransient,'Hcount_plugin_info_TTonic',Hcount_plugin_TTonic,'Hcount_plugin_info_TOff',Hcount_plugin_TOff,'max_info_tpmc_idx_TTransient',max_info_tpmc_idx_TTransient,'max_info_tpmc_TTransient',max_info_tpmc_TTransient,'max_info_tpmc_idx_TTonic',max_info_tpmc_idx_TTonic,'max_info_tpmc_TTonic',max_info_tpmc_TTonic,'max_info_tpmc_idx_TOff',max_info_tpmc_idx_TOff,'max_info_tpmc_TOff',max_info_tpmc_TOff,'Hcount_tpmc_info_TTransient',Hcount_plugin_TTransient,'Hcount_tpmc_info_TTonic',Hcount_tpmc_TTonic,'Hcount_tpmc_info_TOff',Hcount_tpmc_TOff,'max_info_jack_idx_TTransient',max_info_jack_idx_TTransient,'max_info_jack_TTransient',max_info_jack_TTransient,'Hcount_jack_info_TTransient',Hcount_jack_TTransient,'HBias_TTransient',HBias_TTransient,'HBias_std_TTransient',HBias_std_TTransient,'max_info_jack_idx_TTonic',max_info_jack_idx_TTonic,'max_info_jack_TTonic',max_info_jack_TTonic,'Hcount_jack_info_TTonic',Hcount_jack_TTonic,'HBias_TTonic',HBias_TTonic,'HBias_std_TTonic',HBias_std_TTonic,'max_info_jack_idx_TOff',max_info_jack_idx_TOff,'max_info_jack_TOff',max_info_jack_TOff,'Hcount_jack_info_TOff',Hcount_jack_TOff,'HBias_TOff',HBias_TOff,'HBias_std_TOff',HBias_std_TOff,'info_unjkTransient',info_unjkTransient,'info_unjkTonic',info_unjkTonic,'info_unjkOff',info_unjkOff,'HfitTransient',HfitTransient,'HfitTonic',HfitTonic,'HfitOff',HfitOff,'qTransient',qTransient,'qTonic',qTonic,'qOff',qOff,'tempPrecisionTransient',tempPrecisionTransient,'tempPrecisionTonic',tempPrecisionTonic,'tempPrecisionOff',tempPrecisionOff,'percTempTransient',percTempTransient,'percTempTonic',percTempTonic,'percTempOff',percTempOff);

            metric_analysis = struct('XTransient', XTransient, 'XTonic', XTonic, 'XOff', XOff, 'out_TTransient', out_TTransient, 'out_TTonic', out_TTonic, 'out_TOff', out_TOff, 'max_info', max_info, 'opts', opts, 'Y_TTransient', Y_TTransient, 'Y_TTonic', Y_TTonic, 'Y_TOff', Y_TOff, 'SHUF_TTransient', SHUF_TTransient, 'SHUF_TTonic', SHUF_TTonic,'SHUF_TOff', SHUF_TOff);

            save(strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-metric_analysis-graph-5'),'metric_analysis');  

            disp('PLOTTING INFO TIMING-INTERVAL');
            % INFO TIMING-INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f(A) = figure;
            %plot(1:length(opts.shift_cost),info_tpmc_TTransient,'bo');
            errorbar(1:length(opts.shift_cost),info_tpmc_TTransient,2*info_jk_semTransient,'bo');
            hold on;
            
            %plot(1:length(opts.shift_cost),info_tpmc_TTonic,'ro');
            errorbar(1:length(opts.shift_cost),info_tpmc_TTonic,2*info_jk_semTonic,'ro');
            
            %plot(1:length(opts.shift_cost),info_tpmc_TOff,'go');
            errorbar(1:length(opts.shift_cost),info_tpmc_TOff,2*info_jk_semOff,'go');
            hold off;

            set(gca,'xtick',1:length(opts.shift_cost));
            set(gca,'xticklabel',opts.shift_cost);
            set(gca,'xlim',[1 length(opts.shift_cost)]);
            set(gca,'ylim',[-0.5 2.5]);

            xlabel('Temporal precision (1/sec)');
            ylabel('Information (bits)');

            legend('Transient','Tonic','Off');

            print(f(A),'-depsc',strcat(filepath,char(registros{c}),'-contrast-',int2str(A),'-plot-info-timing-transient'));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end     
        
    end


    function H = Hfit(x,q)
        
        H = x(1)*((1 + x(2)*q.^x(3))./(1 + x(4)*q.^x(5)));
        
    end

toc

end

