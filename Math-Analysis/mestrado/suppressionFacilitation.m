
function suppressionFacilitation(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions,CRF,nCRF)

tic

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

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

     for i=1:nConditions
        
       trials_spikes = spike_times(trials_label(i).label,:);
        
       nTrials(i) = size(trials_spikes,1);

     end

    minTrials = min(nTrials);

    for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).label,:);
        
       nTrials = size(trials_spikes,1);
       
       for j=1:minTrials
          
           trial = trials_spikes(j,:);
           
           trial = trial(trial>0);
           
           trial = trial(trial>(start_time/1000) & trial<(end_time/1000));
           
           nBins = (end_time - start_time)/bin_size;
           
           for k=1:nBins
          
                spikes = length(trial(trial>=((k-1)*bin_size/1000 + start_time/1000) & trial<(k*bin_size/1000 + start_time/1000)));
           
                condition(i).allTrials(j).binRateTrial(k) = spikes/(bin_size/1000);

           end
       
           condition(i).meanFRTrial(j) = mean(condition(i).allTrials(j).binRateTrial);
       end
       
       spikes_vector = reshape(trials_spikes.',[],1);

       spikes_vector = spikes_vector.';
       
       spikes_vector = sort(spikes_vector);
       
       spikes_vector = spikes_vector(spikes_vector>0);

       spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));
       
       nBins = (end_time - start_time)/bin_size;
       
       suppressionAnalysis.condition(i).allSpikes = length(spikes_vector);
       
       suppressionAnalysis.condition(i).meanSpikesPerTrial = suppressionAnalysis.condition(i).allSpikes / minTrials;
       
       suppressionAnalysis.condition(i).meanRate = suppressionAnalysis.condition(i).meanSpikesPerTrial / ( (end_time/1000) - (start_time/1000) );
       
       for k=1:nBins
          
           spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
           suppressionAnalysis.conditions(i).binRate(k) = spikes/(nTrials*bin_size/1000);

       end
       
    end
    
    hFRCRF = lillietest(condition(CRF).meanFRTrial);
    hFRNonCRF = lillietest(condition(nCRF).meanFRTrial);
    if nConditions == 3, hFRextraNonCRF = lillietest(condition(3).meanFRTrial); end
        
    [pFRNon,hFRNon] = signrank(condition(1).meanFRTrial,condition(nCRF).meanFRTrial);
    if nConditions == 3, [pFRExtra,hFRExtra] = signrank(condition(CRF).meanFRTrial,condition(3).meanFRTrial); end
    
    suppressionAnalysis.conditions(CRF).meanFRTrial = condition(CRF).meanFRTrial;
    suppressionAnalysis.conditions(nCRF).meanFRTrial = condition(nCRF).meanFRTrial;
    if nConditions == 3, suppressionAnalysis.conditions(3).meanFRTrial = condition(3).meanFRTrial; end
    
    suppressionAnalysis.conditions(CRF).HFRlillietest = hFRCRF;
    suppressionAnalysis.conditions(nCRF).HFRlillietest = hFRNonCRF;
    if nConditions == 3, suppressionAnalysis.conditions(3).HFRlillietest = hFRextraNonCRF; end
    
    suppressionAnalysis.pFRNon = pFRNon;
    suppressionAnalysis.hFRNon = hFRNon;
    if nConditions == 3
        
        suppressionAnalysis.pFRExtra = pFRExtra;
        suppressionAnalysis.hFRExtra = hFRExtra;
    
    end
    
    hrateCRF = lillietest(suppressionAnalysis.conditions(CRF).binRate);
    hrateNonCRF = lillietest(suppressionAnalysis.conditions(nCRF).binRate);
    if nConditions == 3, hrateextraNonCRF = lillietest(suppressionAnalysis.conditions(3).binRate); end
        
    [prateNon,hrateNon] = signrank(suppressionAnalysis.conditions(CRF).binRate,suppressionAnalysis.conditions(2).binRate);
    if nConditions == 3, [prateExtra,hrateExtra] = signrank(suppressionAnalysis.conditions(CRF).binRate,suppressionAnalysis.conditions(3).binRate); end
        
    suppressionAnalysis.conditions(CRF).Hratelillietest = hrateCRF;
    suppressionAnalysis.conditions(nCRF).Hratelillietest = hrateNonCRF;
    if nConditions == 3, suppressionAnalysis.conditions(3).Hratelillietest = hrateextraNonCRF; end
    
    suppressionAnalysis.prateNon = prateNon;
    suppressionAnalysis.hrateNon = hrateNon;
    if nConditions == 3
        
        suppressionAnalysis.prateExtra = prateExtra;
        suppressionAnalysis.hrateExtra = hrateExtra;
        
    end
    
    for k=1:nBins
       
        condition(CRF).trialsBins = [];
        condition(nCRF).trialsBins = [];
        if nConditions == 3, condition(3).trialsBins = []; end
        
        for i=1:minTrials
            
            condition(CRF).trialsBins(i) = condition(CRF).allTrials(i).binRateTrial(k);
            condition(nCRF).trialsBins(i) = condition(nCRF).allTrials(i).binRateTrial(k);
            if nConditions == 3, condition(3).trialsBins(i) = condition(3).allTrials(i).binRateTrial(k); end
            
        end
        
        hrateBinCRF(k) = lillietest(condition(CRF).trialsBins);
        hrateBinNonCRF(k) = lillietest(condition(nCRF).trialsBins);
        if nConditions == 3, hrateBinExtraNonCRF(k) = lillietest(condition(3).trialsBins); end
        
        [prateBinNon(k),hrateBinNon(k)] = signrank(condition(1).trialsBins,condition(nCRF).trialsBins);
        if nConditions == 3, [prateBinExtraNon(k),hrateBinExtraNon(k)] = signrank(condition(CRF).trialsBins,condition(3).trialsBins); end
        
        meanRateBinCRF(k) = mean(condition(CRF).trialsBins);
        meanRateBinNon(k) = mean(condition(nCRF).trialsBins);
        if nConditions == 3, meanRateBinExtraNon(k) = mean(condition(3).trialsBins); end

    end
  
    suppressionAnalysis.prateBinNon = prateBinNon;
    suppressionAnalysis.hrateBinNon = hrateBinNon;
    if nConditions == 3 
        
        suppressionAnalysis.prateBinExtraNon = prateBinExtraNon;
        suppressionAnalysis.hrateBinExtraNon = hrateBinExtraNon;
        
    end
    
    suppressionAnalysis.meanRateBinCRF = meanRateBinCRF;
    suppressionAnalysis.meanRateBinNon = meanRateBinNon;
    if nConditions == 3, suppressionAnalysis.meanRateBinExtraNon = meanRateBinExtraNon; end
    
    suppressionAnalysis.hrateBinCRF = hrateBinCRF;
    suppressionAnalysis.hrateBinNonCRF = hrateBinNonCRF;
    if nConditions == 3, suppressionAnalysis.hrateBinExtraNonCRF = hrateBinExtraNonCRF; end
        
    for k=1:nBins
       
        suppressionAnalysis.NonCRFxCRF(k) = suppressionAnalysis.conditions(nCRF).binRate(k) - suppressionAnalysis.conditions(CRF).binRate(k); 
        
        suppressionAnalysis.RateModulationNonCRF(k) = suppressionAnalysis.conditions(nCRF).binRate(k)/suppressionAnalysis.conditions(CRF).binRate(k);
        
        if nConditions == 3
            
            suppressionAnalysis.ExtraNonCRFxCRF(k) = suppressionAnalysis.conditions(3).binRate(k) - suppressionAnalysis.conditions(CRF).binRate(k);
            
            suppressionAnalysis.RateModulationExtraNonCRF(k) = suppressionAnalysis.conditions(3).binRate(k)/suppressionAnalysis.conditions(CRF).binRate(k);
        
        end
        
    end
    
%     indexZeroBinsNonCRF = find(suppressionAnalysis.RateModulationNonCRF(suppressionAnalysis.RateModulationNonCRF == 0));
%     indexINFBinsNonCRF = find(suppressionAnalysis.RateModulationNonCRF(isinf(suppressionAnalysis.RateModulationNonCRF)));
%     indexNaNBinsNonCRF = find(suppressionAnalysis.RateModulationNonCRF(isnan(suppressionAnalysis.RateModulationNonCRF)));
    
    RateModulationNonCRFfinite = suppressionAnalysis.RateModulationNonCRF(isfinite(suppressionAnalysis.RateModulationNonCRF));
    RateModulationNonCRFfiniteF = RateModulationNonCRFfinite(isfinite(RateModulationNonCRFfinite));
    
    suppressionAnalysis.meanRateModulationNonCRF = mean(RateModulationNonCRFfinite);
    suppressionAnalysis.conditions(CRF).meanFiringRate = mean(suppressionAnalysis.conditions(CRF).binRate);
    suppressionAnalysis.conditions(nCRF).meanFiringRate = mean(suppressionAnalysis.conditions(nCRF).binRate);
    
    if nConditions == 3
        
        suppressionAnalysis.conditions(3).meanFiringRate = mean(suppressionAnalysis.conditions(3).binRate);
    
%         indexZeroBinsExtraNonCRF = find(suppressionAnalysis.RateModulationExtraNonCRF(suppressionAnalysis.RateModulationExtraNonCRF == 0));
%         indexINFBinsExtraNonCRF = find(suppressionAnalysis.RateModulationExtraNonCRF(isinf(suppressionAnalysis.RateModulationExtraNonCRF)));
%         indexNaNBinsExtraNonCRF = find(suppressionAnalysis.RateModulationExtraNonCRF(isnan(suppressionAnalysis.RateModulationExtraNonCRF)));
    
        RateModulationExtraNonCRFfinite = suppressionAnalysis.RateModulationExtraNonCRF(isfinite(suppressionAnalysis.RateModulationExtraNonCRF));
        RateModulationExtraNonCRFfinite = RateModulationExtraNonCRFfinite(isfinite(RateModulationExtraNonCRFfinite));
       
        suppressionAnalysis.meanRateModulationExtraNonCRF = mean(RateModulationExtraNonCRFfinite);
    
    end
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    mkdir(filepath,'suppression');
    filepath = strcat(filepath,'/','suppression','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-suppression');
       
    save(filepath,'suppressionAnalysis');
    
    
        f1 = figure;

        plot(1:nBins,suppressionAnalysis.NonCRFxCRF)
        ylabel(strcat('Suppression x Facilitation: NonCRFxCRF '));
        xlabel('Time');
        %set (gcf, 'WindowButtonMotionFcn', {@mouseMove, video_index, nBins});
        print(f1,'-djpeg',strcat(filepath,'-SuppressionxFacilitationNonCRFxCRF.bmp'));
        
         
        f2 = figure;
        
        plot(1:size(suppressionAnalysis.RateModulationNonCRF,2),suppressionAnalysis.RateModulationNonCRF)
        hold on;
        
%         hrateBinCRF(indexZeroBinsNonCRF) = [];
%         hrateBinNonCRF(indexZeroBinsNonCRF) = [];
%         prateBinNon(indexZeroBinsNonCRF) = [];
%         hrateBinNon(indexZeroBinsNonCRF) = [];
%         
%         hrateBinCRF(indexINFBinsNonCRF) = [];
%         hrateBinNonCRF(indexINFBinsNonCRF) = [];
%         prateBinNon(indexINFBinsNonCRF) = [];
%         hrateBinNon(indexINFBinsNonCRF) = [];
%         
%         hrateBinCRF(indexNaNBinsNonCRF) = [];
%         hrateBinNonCRF(indexNaNBinsNonCRF) = [];
%         prateBinNon(indexNaNBinsNonCRF) = [];
%         hrateBinNon(indexNaNBinsNonCRF) = [];
        
        for i=1:size(suppressionAnalysis.RateModulationNonCRF,2)
            
          if ((hrateBinCRF(i) == 1) && (hrateBinNonCRF(i) == 1) && (hrateBinNon(i) == 1) && (prateBinNon(i) < 0.05))
        
              if suppressionAnalysis.RateModulationNonCRF(i) < 1;

                left = i;
                right = i + 0.5;
                top = suppressionAnalysis.RateModulationNonCRF(i);
                bottom = 0;
                x = [left left right right];
                y = [bottom top top bottom];
                  
                  fill(x,y, 'b','EdgeColor','none');
                  
              else
                  
                left = i;
                right = i + 0.5;
                top = suppressionAnalysis.RateModulationNonCRF(i);
                bottom = 1;
                x = [left left right right];
                y = [bottom top top bottom];

                  fill(x,y, 'g','EdgeColor','none');
                  
              end
              
          end
          
        end     

            p1 = [1 1];
            p2 = [0 size(suppressionAnalysis.RateModulationNonCRF,2)];
            plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);
        
        ylabel('{Taxa de Modula\c{c}\~ao}', 'interpreter', 'latex','FontSize',30);
        xlabel('Bins','FontSize',50);
        %title('{Modula\c{c}\~ao do Contorno 3x o Tamanho do Campo}', 'interpreter', 'latex','FontSize',30)
        %set (gcf, 'WindowButtonMotionFcn', {@mouseMove, video_index, nBins});
        print(f2,'-djpeg',strcat(filepath,'-RateModulationNonCRF.bmp'));
     
         [uniques,numUniques] = count_unique(suppressionAnalysis.RateModulationNonCRF);
         
         uniques = uniques(isfinite(uniques));
   
         uniques = uniques(uniques>0);
         
%          for i=1:size(uniques,1)
% 
%               allUniquesNonCRF(round(uniques(i)*100)) = numUniques(i);
%              
%          end

        [uniques, sortIdx] = sort(uniques);
        numUniques = numUniques(sortIdx);

        frequencias = numUniques;
        for i=2:length(frequencias)

            for k=1:i-1

                frequencias(i) = frequencias(i) + frequencias(k);

            end

        end

        frequencias = frequencias./max(frequencias);

         f3 = figure;
         
         plot(uniques,frequencias);
         %xlim([0 1400]);
         xlabel('{N\''iveis de Taxa de Modula\c{c}\~ao}','interpreter','latex','FontSize',30);
         ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
         title('{Distribui\c{c}\~ao da Taxa de Modula\c{c}\~ao}','interpreter','latex','FontSize',30);
         print(f3,'-djpeg',strcat(filepath,'-RateModulationHistDistributionNonCRF.bmp'));

        if nConditions == 3
            
            
            f4 = figure;

            plot(1:nBins,suppressionAnalysis.ExtraNonCRFxCRF)
            ylabel(strcat('Suppression x Facilitation: ExtraNonCRFxCRF '));
            xlabel('Time');
            %set (gcf, 'WindowButtonMotionFcn', {@mouseMove, video_index, nBins});
            print(f4,'-djpeg',strcat(filepath,'-SuppressionxFacilitationExtraNonCRFxCRF.bmp'));
            
            f5 = figure;
        
            plot(1:size(suppressionAnalysis.RateModulationExtraNonCRF,2),suppressionAnalysis.RateModulationExtraNonCRF)
            hold on;
            
        hrateBinCRF = suppressionAnalysis.hrateBinCRF;

%         hrateBinCRF(indexZeroBinsExtraNonCRF) = [];
%         hrateBinExtraNonCRF(indexZeroBinsExtraNonCRF) = [];
%         prateBinExtraNon(indexZeroBinsExtraNonCRF) = [];
%         hrateBinExtraNon(indexZeroBinsExtraNonCRF) = [];
%         
%         hrateBinCRF(indexINFBinsExtraNonCRF) = [];
%         hrateBinExtraNonCRF(indexINFBinsExtraNonCRF) = [];
%         prateBinExtraNon(indexINFBinsExtraNonCRF) = [];
%         hrateBinExtraNon(indexINFBinsExtraNonCRF) = [];
%         
%         hrateBinCRF(indexNaNBinsExtraNonCRF) = [];
%         hrateBinExtraNonCRF(indexNaNBinsExtraNonCRF) = [];
%         prateBinExtraNon(indexNaNBinsExtraNonCRF) = [];
%         hrateBinExtraNon(indexNaNBinsExtraNonCRF) = [];

        for i=1:size(suppressionAnalysis.RateModulationExtraNonCRF,2)
            
          if ((hrateBinCRF(i) == 1) && (hrateBinExtraNonCRF(i) == 1) && (hrateBinExtraNon(i) == 1) && (prateBinExtraNon(i) < 0.05))
        
              if suppressionAnalysis.RateModulationExtraNonCRF(i) < 1;

                left = i;
                right = i + 0.5;
                top = suppressionAnalysis.RateModulationExtraNonCRF(i);
                bottom = 0;
                x = [left left right right];
                y = [bottom top top bottom];
                  
                  fill(x,y, 'b','EdgeColor','none');
                  
              else
                  
                left = i;
                right = i + 0.5;
                top = suppressionAnalysis.RateModulationExtraNonCRF(i);
                bottom = 1;
                x = [left left right right];
                y = [bottom top top bottom];

                  fill(x,y, 'g','EdgeColor','none');
                  
              end
              
          end
          
        end     
        
            p1 = [1 1];
            p2 = [0 size(suppressionAnalysis.RateModulationExtraNonCRF,2)];
            plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);

            ylabel('{Taxa de Modula\c{c}\~ao}', 'interpreter', 'latex','FontSize',30);
            xlabel('Bins','FontSize',50);
            %title('{Modula\c{c}\~ao do Contorno 6x o Tamanho do Campo}', 'interpreter', 'latex','FontSize',30)
            %set (gcf, 'WindowButtonMotionFcn', {@mouseMove, video_index, nBins});
            print(f5,'-djpeg',strcat(filepath,'-RateModulationExtraNonCRF.bmp'));
        
    

         [uniques,numUniques] = count_unique(suppressionAnalysis.RateModulationExtraNonCRF);

         uniques = uniques(isfinite(uniques));
         
         uniques = uniques(uniques>0);
   
%         for i=1:size(uniques,1)
% 
%             allUniquesExtraNonCRF(round(uniques(i)*100)) = numUniques(i);
%         
%         end

        f6 = figure;

        [uniques, sortIdx] = sort(uniques);
        numUniques = numUniques(sortIdx);

        frequencias = numUniques;
        for i=2:length(frequencias)

            for k=1:i-1

                frequencias(i) = frequencias(i) + frequencias(k);

            end

        end

        frequencias = frequencias./max(frequencias);

        plot(uniques,frequencias);
        %xlim([0 1400]);
        xlabel('{N\''iveis de Taxa de Modula\c{c}\~ao}','interpreter','latex','FontSize',30);
        ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
        title('{Distribui\c{c}\~ao da Taxa de Modula\c{c}\~ao}','interpreter','latex','FontSize',30);
        print(f6,'-djpeg',strcat(filepath,'-RateModulationHistDistributionExtraNonCRF.bmp'));

        end
    
    
toc
    

% function C = mouseMove (object, eventdata, video_index, nBins)
%     
%     clear h;
%     
%     C = get (gca, 'CurrentPoint');
%     title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
%     
%     if ( C(1,1) > 0 ) && ( C(1,1) < nBins )
%         
%           img = imread(strcat('v',int2str(video_index),'-radius-',int2str(22),'-frame-',int2str(round(C(1,1))),'.jpg'));
%               
%           %newimg = imresize(img,[12 12]);
%           
%           hold on;
%           
%           h = image(C(1,1),40,img);
%           
%     end
%     
% end


end
       