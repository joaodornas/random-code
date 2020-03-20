
function binOptimalPSTHkernel_Ana(allConditions,AtividadeEspontanea) %binOptimalPSTHkernel(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

tic
    
A = exist('PSTHKernel','dir');

if A ~=7, mkdir('PSTHKernel'); end
    
fileattrib('PSTHKernel','+w');

cells = {'gas004a01_3a' 'gas004a02_1c' 'gas005a01_1b' 'gas005c01_1b' 'gas006a02_1b' 'gas008a01_1c' 'gas009a02_2b' 'gas009a02_2c' 'gas009a04_3a' 'gas010a02_2b' 'gas011b02_1b' 'gas012a02_1b' 'gas013a02_2c' 'gas016a02_3a' 'gas020b02_1b' 'gas024a03_2B' 'gas024a03_2C' 'gas025a02_1B' 'gas025a02_1C' 'gas026a03_3B' 'gas031a02_1B' 'gas034a02_1B' 'gas036a03_2B' 'gas036a03_2C' 'gas037a03_1B' 'gas038b02_1B' 'gas039a02_1C' 'gas040a01_2B' 'gas043b02_2B' 'gas046b02_1B' 'gas049a02_2C' 'gas051a02_1B' 'gas051d02_1B' 'gas053b02_1A' 'gas053b03_2B' 'gas054b01_2B' 'gas055a02_1B' 'gas056a02_1B' 'gas057a02_1B' 'gas057a02_1C' 'gas057b03_2B' 'gas058b02_1A' 'gas059a02_1B' 'gas061a03_2B' 'gas061b02_2B' 'gas061b03_1A' 'gas061b03_1C' 'gas062b02_1B' 'gas063a02_1B' 'gas066b01_1B' 'gas068b02_1A' 'gas071a02_2B' 'gas071a03_2A' 'gas072a05_1B' 'gas075a02_3A' 'gas075b02_2B' 'gas076a02_1B' 'gas081b02_1B' 'gas083a02_1B' 'gas083a04_3B' 'gas084a04_3b' 'gas087a02_3B' 'gas091a02_3a' 'gas091b02_3b' 'gas093a02_2B' 'gas094a02_1B' 'gas095a02_1B' 'gas096a02_1B' 'gas109a03_2b' 'gas110a03_2a'};
stim_start_time = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000]; 
stim_end_time = [5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000];
duration = [7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 7000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000]; 
startCondition = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 16 16 16 16 18 18 18 18 1 18 1 1 18 18 18 18 1 1];
endCondition = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 30 30 30 30 34 34 34 34 10 34 17 17 34 34 34 34 16 16];
totalConditions = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 30 30 30 30 34 34 34 34 10 34 17 17 34 34 34 34 32 32];

for c=1:length(cells)
    
    disp('BEGIN');
    disp(strcat('Loading cell-',int2str(c)));
    disp(strcat('Cell name:',char(cells(c))));
    
    if allConditions == 0
        
        firstCondition = startCondition(c);
        lastCondition = endCondition(c);
        
    elseif allConditions == 1
        
        firstCondition = 1;
        lastCondition = totalConditions(c);
        
    end
    
    if AtividadeEspontanea == 0
        
        start_time = stim_start_time(c);
        end_time = stim_end_time(c);
        
    elseif AtividadeEspontanea == 1
        
        start_time = 0;
        end_time = duration(c);
        
    end
    
    %%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Load Spass...');

        Spass = load(strcat(char(cells(c)),'.mat'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Load Conditions Trials Labels...');

        for i=firstCondition:lastCondition

            trials_label(i).labels = find(Spass.stimIds == i); 

        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Set Resolution...');

        spike_times = Spass.spike_times;

        spike_times = spike_times./ 32000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        histData = struct('trials_spikes',[],'spikes_vector',[],'psth',struct('Cost',[],'NBinsTested',[],'binHeight',[],'binSize',0,'binHeightNormal',[]),'kernel',struct('optimalBinWidth',0,'Cost',[],'WBinsTested',[], 'density',[], 'timePoints', []),'nTrials',0); 

    %%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Begin Conditions Loop...');

        for i=firstCondition:lastCondition

           disp(strcat('Begin Condition . ',int2str(i))); 

           trials_spikes = spike_times(trials_label(i).labels,:);

           nTrials = size(trials_spikes,1);

           spikes_vector = reshape(trials_spikes.',[],1);

           spikes_vector = spikes_vector.';

           spikes_vector = sort(spikes_vector);

           spikes_vector = spikes_vector(spikes_vector>0);

           spikes_vector = spikes_vector(spikes_vector>=(start_time/1000) & spikes_vector<=(end_time/1000));

           NbinMax = round((end_time/1000 - start_time/1000) * 1000/2);

           bandwidths = (2/1000):(2/1000):((end_time - start_time)/1000/2);

           %%%  PSTH COST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...PSTH Cost - Shinomoto');

           [histData.condition(i).psth.optN, histData.condition(i).psth.Cost, histData.condition(i).psth.NBinsTested] = sshist(spikes_vector,2:NbinMax);

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %%%  KERNEL COST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...Kernel Cost - Shinomoto'); 

           time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * 1000);

           [histData.condition(i).kernel.density, histData.condition(i).kernel.timePoints, histData.condition(i).kernel.optimalBinWidth, histData.condition(i).kernel.WBinsTested, histData.condition(i).kernel.Cost, histData.condition(i).kernel.confb95] = ssvkernel(spikes_vector,time_points);

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %%%  SAVE VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...save variables');

           histData.condition(i).trials_spikes = trials_spikes;
           histData.condition(i).spikes_vector = sort(spikes_vector);
           histData.condition(i).psth.NbinMax = NbinMax;
           histData.condition(i).nTrials = nTrials;

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %%%  MAKE HISTOGRAM DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...make histogram data');

           binSize = (end_time/1000 - start_time/1000)/histData.condition(i).psth.optN;

           for k=1:histData.condition(i).psth.optN

               spikes = length(spikes_vector(spikes_vector>=((k-1)*binSize + start_time/1000) & spikes_vector<(k*binSize + start_time/1000)));

               binHeight(k) = spikes/(nTrials*binSize);

           end

           %%%  SAVE AND CLEAR VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...save and clear variables');

           histData.condition(i).psth.binHeight = binHeight;
           histData.condition(i).psth.binSize = binSize;

           histData.condition(i).psth.binHeightNormal = histData.condition(i).psth.binHeight./max(histData.condition(i).psth.binHeight);

           clear trials_spikes
           clear spikes_vector
           clear nTrials
           clear NbinMax bandwidths
           clear binHeight binSize

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           disp(strcat('End Condition . ',int2str(i)));

            histData.cell_label = strcat(char(cells(c)));

            histData.condition(i).timePoints = length(histData.condition(i).kernel.timePoints);
            histData.condition(i).timeHistoryBegin = start_time/1000 + 20/1000 ;
            histData.condition(i).timeHistoryEnd = start_time/1000 + 150/1000 + 20/1000;
            histData.condition(i).totalTime = end_time/1000 - start_time/1000 ;
            histData.condition(i).timeRatio = histData.condition(i).timePoints / histData.condition(i).totalTime ;

            histData.condition(i).Ypico = max(histData.condition(i).kernel.density(histData.condition(i).timeHistoryBegin*histData.condition(i).timeRatio:histData.condition(i).timeHistoryEnd*histData.condition(i).timeRatio));

            for k=histData.condition(i).timeHistoryBegin*histData.condition(i).timeRatio:histData.condition(i).timeHistoryEnd*histData.condition(i).timeRatio

                if histData.condition(i).kernel.density(k) == histData.condition(i).Ypico

                    histData.condition(i).Xpico = k;

                end

            end

            histData.condition(i).ae = median(histData.condition(i).kernel.density(1:histData.condition(i).timeHistoryBegin*histData.condition(i).timeRatio));

            histData.condition(i).A = ( (histData.condition(i).Ypico - histData.condition(i).ae) / 2 ) + histData.condition(i).ae ;
            histData.condition(i).kernel.density = histData.condition(i).kernel.density';
            histData.condition(i).V = histData.condition(i).kernel.density(histData.condition(i).timeHistoryBegin*histData.condition(i).timeRatio:histData.condition(i).Xpico);
            histData.condition(i).XM = dsearchn(histData.condition(i).V,histData.condition(i).A);
            histData.condition(i).XM = histData.condition(i).XM + histData.condition(i).timeHistoryBegin*histData.condition(i).timeRatio;
            histData.condition(i).latencyPoints = histData.condition(i).kernel.timePoints(histData.condition(i).XM);
           
        end


    disp('End Conditions Loop');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%% PLOT DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Begin Plot Data Loop...');    

        for i=firstCondition:lastCondition

           disp(strcat('Begin Condition . ',int2str(i)));

           %%%   KERNEL DENSITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...Kernel density');

           %[histData(i).kernel.KY,histData(i).kernel.KX,histData(i).kernel.optimalBinWidth] = sskernel(histData(i).spikes_vector,t,histData(i).kernel.optimalBinWidth);



           more = length(histData.condition(i).kernel.density(histData.condition(i).kernel.density>end_time/1000));
           less = length(histData.condition(i).kernel.timePoints(histData.condition(i).kernel.timePoints<start_time/1000));
           histData.condition(i).kernel.density = histData.condition(i).kernel.density(less+1:end-more);
           histData.condition(i).kernel.timePoints = histData.condition(i).kernel.timePoints(less+1:end-more);

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %        %%%   ERROR DENSITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        disp('...error density');
    %        
    %        Gauss = @(s,w) 1/sqrt(2*pi)/w*exp(-s.^2/2/w^2);
    %        nPoints = length(histData.condition(i).kernel.density);
    %        
    %        for k=1:nPoints
    %        
    %            histData.condition(i).kernel.error(k) = std( Gauss(histData.condition(i).spikes_vector-histData.condition(i).kernel.density(k),histData.condition(i).kernel.optimalBinWidth) );
    %                       
    %        end
    %        
    %        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %%%   PLOT FIGURES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           disp('...plot figures');

           %%%   FIGURE: ALL RASTER + KERNEL    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           f(i) = figure;
           disp(strcat('figure:',int2str(i),'.1'));

           disp('...all raster plot');
           subplot(2,1,1);

           for k=1:size(histData.condition(i).trials_spikes,1)

               spikes = histData.condition(i).trials_spikes(k,:);
               spikes = spikes(spikes>0);
               spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
               spikes = sort(spikes);

               for j=1:size(spikes,2)

                   plot([spikes(j) spikes(j)],[(k-0.4) (k)],'b');
                   hold on;

               end

           end

           disp('...kernel density function');
           subplot(2,1,2);

           plot(histData.condition(i).kernel.timePoints,histData.condition(i).kernel.density);
           title(strcat(char(cells(c)),'-condition-',int2str(i)));

           %%%   FIGURE: ALL RASTER + PSTH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           h(i) = figure;
           disp(strcat('figure:',int2str(i),'.3'));

           disp('...all raster plot');
           subplot(2,1,1);

           for k=1:size(histData.condition(i).trials_spikes,1)

               spikes = histData.condition(i).trials_spikes(k,:);
               spikes = spikes(spikes>0);
               spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
               spikes = sort(spikes);

               for j=1:size(spikes,2)

                   plot([spikes(j) spikes(j)],[(k-0.4) (k)],'b');
                   hold on;

               end

           end

           disp('...psth');
           subplot(2,1,2);

           bar(1:histData.condition(i).psth.optN,histData.condition(i).psth.binHeightNormal);
           title(strcat(char(cells(c)),'-condition-',int2str(i)));

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           disp('...saving pictures');
           
%            A = exist(strcat('/','PSTHKernel','/',char(cells(c))),'dir');
%            
%            if A ~=7, mkdir('PSTHKernel',char(cells(c))); end
%            
%            fileattrib(strcat('PSTHKernel','/',char(cells(c))),'+w');
%     
%            filepath = strcat('/','PSTHKernel','/',strcat(char(cells(c))),'/');
                    
             filepath = strcat('/','PSTHKernel','/');
                    
%            
%            A = exist(strcat('/','PSTHKernel','/',char(cells(c)),'/','kernel'),'dir');
%            
%            if A ~=7, mkdir(strcat('PSTHKernel/',char(cells(c))),'kernel'); end
%            
%            fileattrib(strcat('PSTHKernel','/',char(cells(c)),'/','kernel'),'+w');
%            
%            A = exist(strcat('/','PSTHKernel','/',char(cells(c)),'/','psth'),'dir');
%            
%            if A ~=7, mkdir(strcat('PSTHKernel/',char(cells(c))),'psth'); end
%            
%            fileattrib(strcat('PSTHKernel','/',char(cells(c)),'/','psth'),'+w');
%            
%            filepathPSTH = strcat(filepath,'psth','/');
%            filepathKERNEL = strcat(filepath,'kernel','/');
         
           namePSTH = strcat(char(cells(c)),'-raster-all-psth-condition-');
           nameKERNEL = strcat(char(cells(c)),'-raster-all-kernel-condition-');
           
           computerPath = '/Users/joaodornas/Documents/_Research';
           
           print(f(i),'-depsc',strcat(computerPath,filepath,nameKERNEL,int2str(i)));
           print(h(i),'-depsc',strcat(computerPath,filepath,namePSTH,int2str(i)));

           disp(strcat('End Condition . ',int2str(i)));

           close all;
           
        end

    disp('End Plot Data Loop');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Save data...');
        
    save(strcat(computerPath,filepath,char(cells(c))),'histData');
    
    clear histData;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('END');
    
end
    
toc
    
end