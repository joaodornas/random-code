
function binOptimalPSTHkernel(date,site_index,channel,registro,video_index,start_time,end_time,delay,bin_size,nConditions)

tic

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    for i=1:nConditions
   
        trials_label(i).labels = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times;
    
    spike_times = spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %histData = struct('trials_spikes',[],'spikes_vector',[],'psth',struct('Cost',[],'NBinsTested',[],'binHeight',[],'binSize',0,'binHeightNormal',[]),'kernel',struct('optimalBinWidth',0,'Cost',[],'WBinsTested',[], 'density',[], 'timePoints', []),'nTrials',0); 
    
%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

for invert_movie=1:2
    
    for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).labels,:);
        
       nTrials = size(trials_spikes,1);
       
       spikes_vector = reshape(trials_spikes.',[],1);

       spikes_vector = spikes_vector.';
       
       spikes_vector = sort(spikes_vector);
       
       spikes_vector = spikes_vector(spikes_vector>0);

       spikes_vector = spikes_vector(spikes_vector>=(start_time/1000) & spikes_vector<=(end_time/1000));
       
       if (i == 2) && (invert_movie == 2)
           
            nBins = (end_time - start_time)/bin_size;
            
            disp('...invert spike train .');
            
            spikes_inverted = invertSpikeTrain(spikes_vector,start_time,nBins,bin_size);
            spikes_vector = [];
            spikes_vector = spikes_inverted.spikes;
            
        end
       
       NbinMax = round((end_time/1000 - start_time/1000) * 1000/2);
       
       bandwidths = (2/1000):(2/1000):((end_time - start_time)/1000/2);
       
       %%%  PSTH COST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       disp('...PSTH Cost - Shinomoto');
       
       if (invert_movie == 1)
           
            [histData(invert_movie).condition(i).psth.optN, histData(invert_movie).condition(i).psth.Cost, histData(invert_movie).condition(i).psth.NBinsTested] = sshist(spikes_vector,2:NbinMax);
       
       else
           
           if (i == 2)
               
              [histData(invert_movie).condition(i).psth.optN, histData(invert_movie).condition(i).psth.Cost, histData(invert_movie).condition(i).psth.NBinsTested] = sshist(spikes_vector,2:NbinMax);
               
           else
               
               histData(2).condition(i) = histData(1).condition(i);
               
           end
           
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%  KERNEL COST  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       disp('...Kernel Cost - Shinomoto'); 
      
       time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * 1000);
       
       if (invert_movie == 1)
           
             [histData(invert_movie).condition(i).kernel.density, histData(invert_movie).condition(i).kernel.timePoints, histData(invert_movie).condition(i).kernel.optimalBinWidth, histData(invert_movie).condition(i).kernel.WBinsTested, histData(invert_movie).condition(i).kernel.Cost, histData(invert_movie).condition(i).kernel.confb95] = ssvkernel(spikes_vector,time_points);
       
       else
           
           if (i == 2)
               
                [histData(invert_movie).condition(i).kernel.density, histData(invert_movie).condition(i).kernel.timePoints, histData(invert_movie).condition(i).kernel.optimalBinWidth, histData(invert_movie).condition(i).kernel.WBinsTested, histData(invert_movie).condition(i).kernel.Cost, histData(invert_movie).condition(i).kernel.confb95] = ssvkernel(spikes_vector,time_points);
       
           else
               
                histData(2).condition(i) = histData(1).condition(i);
               
           end
           
       end
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%  SAVE VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       disp('...save variables');
       
       histData(invert_movie).condition(i).trials_spikes = trials_spikes;
       histData(invert_movie).condition(i).spikes_vector = sort(spikes_vector);
       histData(invert_movie).condition(i).psth.NbinMax = NbinMax;
       histData(invert_movie).condition(i).nTrials = nTrials;
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%  MAKE HISTOGRAM DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       disp('...make histogram data');
       
       binSize = (end_time/1000 - start_time/1000)/histData(invert_movie).condition(i).psth.optN;
       
       for k=1:histData(invert_movie).condition(i).psth.optN
          
           spikes = length(spikes_vector(spikes_vector>=((k-1)*binSize + start_time/1000) & spikes_vector<(k*binSize + start_time/1000)));
           
           binHeight(k) = spikes/(nTrials*binSize);

       end
       
       %%%  SAVE AND CLEAR VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       disp('...save and clear variables');
       
       histData(invert_movie).condition(i).psth.binHeight = binHeight;
       histData(invert_movie).condition(i).psth.binSize = binSize;
         
       histData(invert_movie).condition(i).psth.binHeightNormal = histData(invert_movie).condition(i).psth.binHeight./max(histData(invert_movie).condition(i).psth.binHeight);
       
       clear trials_spikes
       clear spikes_vector
       clear nTrials
       clear NbinMax bandwidths
       clear binHeight binSize
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       disp(strcat('End Condition . ',int2str(i)));
       
    end
    
end

disp('End Conditions Loop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

histData(1).cell_label = strcat('v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);

histData(1).timePoints = length(histData(invert_movie).condition(i).kernel.timePoints);
histData(1).timeHistoryBegin = 500/1000 + 20/1000 ;
histData(1).timeHistoryEnd = 500/1000 + delay/1000 + 20/1000;
histData(1).totalTime = end_time/1000 - start_time/1000 ;
histData(1).timeRatio = histData(1).timePoints / histData(1).totalTime ;

histData(1).Ypico = max(histData(1).condition(1).kernel.density(histData(1).timeHistoryBegin*histData(1).timeRatio:histData(1).timeHistoryEnd*histData(1).timeRatio));

for i=histData(1).timeHistoryBegin*histData(1).timeRatio:histData(1).timeHistoryEnd*histData(1).timeRatio
    
    if histData(1).condition(1).kernel.density(i) == histData(1).Ypico
        
        histData(1).Xpico = i;
        
    end
    
end

histData(1).ae = median(histData(1).condition(1).kernel.density(1:histData(1).timeHistoryBegin*histData(1).timeRatio));

histData(1).A = ( (histData(1).Ypico - histData(1).ae) / 2 ) + histData(1).ae ;
histData(1).condition(1).kernel.density = histData(1).condition(1).kernel.density';
V = histData(1).condition(1).kernel.density(histData(1).timeHistoryBegin*histData(1).timeRatio:histData(1).Xpico);
histData(1).XM = dsearchn(V,histData(1).A);
histData(1).XM = histData(1).XM + histData(1).timeHistoryBegin*histData(1).timeRatio;
histData(1).latencyPoints = histData(1).condition(1).kernel.timePoints(histData(1).XM);
histData(1).latencyPoints = histData(1).latencyPoints - 500/1000;

% %%% PLOT DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Begin Plot Data Loop...');
% 
% for invert_movie=1:2    
% 
%     for i=1:nConditions
%         
%        disp(strcat('Begin Condition . ',int2str(i)));
%        
%        %%%   KERNEL DENSITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        disp('...Kernel density');
%        
%        %[histData(i).kernel.KY,histData(i).kernel.KX,histData(i).kernel.optimalBinWidth] = sskernel(histData(i).spikes_vector,t,histData(i).kernel.optimalBinWidth);
%        
%        
%         
%        more = length(histData(invert_movie).condition(i).kernel.density(histData(invert_movie).condition(i).kernel.density>end_time/1000));
%        less = length(histData(invert_movie).condition(i).kernel.timePoints(histData(invert_movie).condition(i).kernel.timePoints<start_time/1000));
%        histData(invert_movie).condition(i).kernel.density = histData(invert_movie).condition(i).kernel.density(less+1:end-more);
%        histData(invert_movie).condition(i).kernel.timePoints = histData(invert_movie).condition(i).kernel.timePoints(less+1:end-more);
% 
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        
%        %%%   ERROR DENSITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        disp('...error density');
%        
%        Gauss = @(s,w) 1/sqrt(2*pi)/w*exp(-s.^2/2/w^2);
%        nPoints = length(histData(invert_movie).condition(i).kernel.density);
%        
%        for k=1:nPoints
%        
%            histData(invert_movie).condition(i).kernel.error(k) = std( Gauss(histData(invert_movie).condition(i).spikes_vector-histData(invert_movie).condition(i).kernel.density(k),histData(invert_movie).condition(i).kernel.optimalBinWidth) );
%                       
%        end
%        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%   PLOT FIGURES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        disp('...plot figures');
%        
       %%%   FIGURE: ALL RASTER + KERNEL    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        f(i) = figure;
%        disp(strcat('figure:',int2str(i),'.1'));
%        
%        disp('...all raster plot');
%        subplot(2,1,1);
%    
%        for k=1:size(histData(invert_movie).condition(i).trials_spikes,1)
% 
%            spikes = histData(invert_movie).condition(i).trials_spikes(k,:);
%            spikes = spikes(spikes>0);
%            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
%            spikes = sort(spikes);
% 
%            for j=1:size(spikes,2)
% 
%                plot([spikes(j) spikes(j)],[(k-0.4) (k)],'b');
%                hold on;
% 
%            end
% 
%        end
%     
%        disp('...kernel density function');
%        subplot(2,1,2);
% 
%        plot(histData(invert_movie).condition(i).kernel.timePoints,histData(invert_movie).condition(i).kernel.density);
       %hold on;
       %plot(histData(i).kernel.KX,histData(i).kernel.KY+histData(i).kernel.error,'r');
       %plot(histData(i).kernel.KX,histData(i).kernel.KY-histData(i).kernel.error,'r');

%        %%%   FIGURE: CONDENSED RASTER + KERNEL   %%%%%%%%%%%%%%%%%%%%%%%%%%
%        g(i) = figure;
%        disp(strcat('figure:',int2str(i),'.2'));
%        
%        disp('...condensed raster plot');
%        subplot(2,1,1);
% 
%        for j=1:size(histData(invert_movie).condition(i).spikes_vector,2)
% 
%             plot([histData(invert_movie).condition(i).spikes_vector(j) histData(invert_movie).condition(i).spikes_vector(j)],[(1-0.4) (1)],'b');
%             hold on;
% 
%        end
% 
%        disp('...kernel density function');
%        subplot(2,1,2);
% 
%        plot(histData(invert_movie).condition(i).kernel.timePoints,histData(invert_movie).condition(i).kernel.density);
%        %hold on;
%        %plot(histData(i).kernel.KX,histData(i).kernel.KY+histData(i).kernel.error,'r');
%        %plot(histData(i).kernel.KX,histData(i).kernel.KY-histData(i).kernel.error,'r');
       
       %%%   FIGURE: ALL RASTER + PSTH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        h(i) = figure;
%        disp(strcat('figure:',int2str(i),'.3'));
%        
%        disp('...all raster plot');
%        subplot(2,1,1);
%    
%        for k=1:size(histData(invert_movie).condition(i).trials_spikes,1)
% 
%            spikes = histData(invert_movie).condition(i).trials_spikes(k,:);
%            spikes = spikes(spikes>0);
%            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
%            spikes = sort(spikes);
%            
%            for j=1:size(spikes,2)
% 
%                plot([spikes(j) spikes(j)],[(k-0.4) (k)],'b');
%                hold on;
% 
%            end
% 
%        end
%        
%        disp('...psth');
%        subplot(2,1,2);
% 
%        bar(1:histData(invert_movie).condition(i).psth.optN,histData(invert_movie).condition(i).psth.binHeightNormal);

%        %%%   FIGURE: CONDENSED RASTER + PSTH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        l(i) = figure;
%        disp(strcat('figure:',int2str(i),'.4'));
% 
%        disp('...condensed raster plot');
%        subplot(2,1,1);
% 
%        for j=1:size(histData(invert_movie).condition(i).spikes_vector,2)
% 
%             plot([histData(invert_movie).condition(i).spikes_vector(j) histData(invert_movie).condition(i).spikes_vector(j)],[(1-0.4) (1)],'b');
%             hold on;
% 
%        end
% 
%        disp('...psth');
%        subplot(2,1,2);
% 
%        bar(1:histData(invert_movie).condition(i).psth.optN,histData(invert_movie).condition(i).psth.binHeightNormal);
   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %disp('...saving pictures');
       
       %print(f(i),'-dbmp',strcat(filepath,nameKERNEL,int2str(i),'.jpg'));
       %print(g(i),'-dbmp',strcat(filepath,'-raster-condensed-kernel-condition-',int2str(i),'.bmp'));
       %print(h(i),'-dbmp',strcat(filepath,namePSTH,int2str(i),'.jpg'));
       %print(l(i),'-dbmp',strcat(filepath,'-raster-condensed-psth-condition-',int2str(i),'.bmp'));
%         
%        disp(strcat('End Condition . ',int2str(i)));
%        
%     end
%     
% end

       filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
       %filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
       
       filepath = strcat(filepath,'/','psth-kernel','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
       filepath = strcat(filepath,'-psth-kernel');
    
       namePSTH = '-raster-all-psth-condition-';
       nameKERNEL = '-raster-all-kernel-condition-';
       
       if (i == 2) && (invert_movie == 2)
           
           namePSTH = strcat(namePSTH,'invertido-');
           nameKERNEL = strcat(nameKERNEL,'invertido-');
 
       end

       density1 = histData(1).condition(1).kernel.density./max(histData(1).condition(1).kernel.density); 
       maxY1 = max(histData(1).condition(1).psth.binHeight);
       density1 = density1.*maxY1;
       
       density2 = histData(1).condition(2).kernel.density./max(histData(1).condition(2).kernel.density);
       maxY2 = max(histData(1).condition(2).psth.binHeight);
       density2 = density2.*maxY2;
       
       density3 = histData(2).condition(2).kernel.density./max(histData(2).condition(2).kernel.density);
       maxY3 = max(histData(2).condition(2).psth.binHeight);
       density3 = density3.*maxY3;
       
       
       maxY = max([maxY1 maxY2 maxY3]);
       
       
%        confb95 = histData(1).condition(1).kernel.confb95;
%        confb951 = confb95(1,:)./max(max(confb95,[],1));
%        confb952 = confb95(2,:)./max(max(confb95,[],2));
%        confb951 = confb951.*maxY;
%        confb952 = confb952.*maxY;
       
       w(1) = figure;
       plot(histData(1).condition(1).kernel.timePoints,density1,'b');
%        hold on;
%        h = fill([histData(1).condition(1).kernel.timePoints,histData(1).condition(1).kernel.timePoints(end:-1:1)],[confb952,confb951(end:-1:1)],'r-.');
%        set(h,'facealpha',.1);
%        hold on;
%        plot(histData(1).condition(1).kernel.timePoints,density1);
       ylim([0 maxY]);
       ylabel('Taxa de Disparo por Segundo');
       xlabel('Tempo (s)');
       print(w(1),'-depsc',strcat(filepath,'-condition-1'));
       

%        confb95 = histData(1).condition(2).kernel.confb95;
%        confb951 = confb95(1,:)./max(max(confb95,[],1));
%        confb952 = confb95(2,:)./max(max(confb95,[],2));
%        confb951 = confb951.*maxY;
%        confb952 = confb952.*maxY;
       
       w(2) = figure;
       plot(histData(1).condition(2).kernel.timePoints,density2,'r');
%        hold on;
%        h = fill([histData(1).condition(2).kernel.timePoints,histData(1).condition(2).kernel.timePoints(end:-1:1)],[confb952,confb951(end:-1:1)],'r-.');
%        set(h,'facealpha',.1);
%        hold on;
%        plot(histData(1).condition(2).kernel.timePoints,density2);
       ylim([0 maxY]);
       ylabel('Taxa de Disparo por Segundo');
       xlabel('Tempo (s)');
       print(w(2),'-depsc',strcat(filepath,'-condition-2-reverso'));
       
%        confb95 = histData(2).condition(2).kernel.confb95;
%        confb951 = confb95(1,:)./max(max(confb95,[],1));
%        confb952 = confb95(2,:)./max(max(confb95,[],2));
%        confb951 = confb951.*maxY;
%        confb952 = confb952.*maxY;
              
       w(3) = figure;
       plot(histData(2).condition(2).kernel.timePoints,density3,'r');
%        hold on;
%        h = fill([histData(2).condition(2).kernel.timePoints,histData(2).condition(2).kernel.timePoints(end:-1:1)],[confb952,confb951(end:-1:1)],'r-.');
%        set(h,'facealpha',.1);
%        hold on;
%        plot(histData(2).condition(2).kernel.timePoints,density3);
       ylim([0 maxY]);
       ylabel('Taxa de Disparo por Segundo');
       xlabel('Tempo (s)');
       print(w(3),'-depsc',strcat(filepath,'-condition-2-invertido'));
       
       
       
       
       density = histData(1).condition(3).kernel.density./max(histData(1).condition(3).kernel.density);
       maxY = max(histData(1).condition(3).psth.binHeight);
       density = density.*maxY;
       
       w(4) = figure;
       plot(histData(1).condition(3).kernel.timePoints,density,'g');
       
       ylim([0 maxY]);
       ylabel('Taxa de Disparo por Segundo');
       xlabel('Tempo (s)');
       print(w(4),'-depsc',strcat(filepath,'-condition-3-ae'));
  
    
% disp('End Plot Data Loop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
    
    save(filepath,'histData');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');
    
toc
    
end