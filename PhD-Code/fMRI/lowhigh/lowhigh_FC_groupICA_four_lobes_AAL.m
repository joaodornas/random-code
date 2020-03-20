function lowhigh_FC_groupICA_four_lobes_AAL


%settings_jan_0805;
settings_elena_2905;

preprocessed_step = 'residual';
%preprocessed_step = 'filtered';

doTheMath(settings,preprocessed_step);

%plotResults(settings,preprocessed_step);

%plotCumulative(settings,preprocessed_step);

%plotSixCumulativeTwoConditions(settings,preprocessed_step);

%plotVarianceHistogram(settings,preprocessed_step);

%plotTimeCourses(settings,preprocessed_step);

end

function doTheMath(settings,preprocessed_step)

ICA_preprocessed = 'groupICA';

nTR = 300;

nROI = 90;

area_lobe_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area_lobe');
area_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area');

for run=1:2
    
    disp(strcat('run:',int2str(run)));
    
    occipital_ic = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','occipital_lobe','\','melodic_mix'));
    parietal_ic = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','parietal_lobe','\','melodic_mix'));
    temporal_ic = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','temporal_lobe','\','melodic_mix'));
    frontal_ic = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','frontal_lobe','\','melodic_mix'));
    
    occipital_var = textread(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','occipital_lobe','\','eigenvalues_percent'),'%s');
    parietal_var = textread(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','parietal_lobe','\','eigenvalues_percent'),'%s');
    temporal_var = textread(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','temporal_lobe','\','eigenvalues_percent'),'%s');
    frontal_var = textread(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','frontal_lobe','\','eigenvalues_percent'),'%s');
    
    ICA_n(run).occipital_var = str2double(occipital_var);
    ICA_n(run).parietal_var = str2double(parietal_var);
    ICA_n(run).temporal_var = str2double(temporal_var);
    ICA_n(run).frontal_var = str2double(frontal_var);
    
    ICA_n(run).occipital_per = [0, ICA_n(run).occipital_var'];
    ICA_n(run).occipital_per = diff(ICA_n(run).occipital_per)*100;
    ICA_n(run).parietal_per = [0, ICA_n(run).parietal_var'];
    ICA_n(run).parietal_per = diff(ICA_n(run).parietal_per)*100;
    ICA_n(run).temporal_per = [0, ICA_n(run).temporal_var'];
    ICA_n(run).temporal_per = diff(ICA_n(run).temporal_per)*100;
    ICA_n(run).frontal_per = [0, ICA_n(run).frontal_var'];
    ICA_n(run).frontal_per = diff(ICA_n(run).frontal_per)*100;
    
    ICA_n(run).nICOccipital = size(occipital_ic,2);
    ICA_n(run).nICParietal = size(parietal_ic,2);
    ICA_n(run).nICTemporal = size(temporal_ic,2);
    ICA_n(run).nICFrontal = size(frontal_ic,2);
    
    high_occipital_ic = occipital_ic(1:nTR,:);
    high_parietal_ic = parietal_ic(1:nTR,:);
    high_temporal_ic = temporal_ic(1:nTR,:);
    high_frontal_ic = frontal_ic(1:nTR,:);
    
    rest_occipital_ic = occipital_ic(nTR+1:end,:);
    rest_parietal_ic = parietal_ic(nTR+1:end,:);
    rest_temporal_ic = temporal_ic(nTR+1:end,:);
    rest_frontal_ic = frontal_ic(nTR+1:end,:);
    
%     allComponents_high = [high_occipital_ic,high_parietal_ic,high_temporal_ic,high_frontal_ic];
%     allComponents_rest = [rest_occipital_ic,rest_parietal_ic,rest_temporal_ic,rest_frontal_ic];
%     
%     [ICA_FC(run).rho_high,ICA_FC(run).pval_high] = corr(allComponents_high);
%     [ICA_FC(run).rho_rest,ICA_FC(run).pval_rest] = corr(allComponents_rest);

    for iROI=1:nROI
        
       disp(strcat('ROI:',int2str(iROI)));
       
       if iROI ~= 87
           
           area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(run),'-','aal2std_',int2str(iROI),'\','melodic_mix'));
           area_var = textread(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(run),'-','aal2std_',int2str(iROI),'\','eigenvalues_percent'),'%s');
           
           ICA_n(run).area_var(iROI).var = str2double(area_var);
           
           ICA_n(run).area_per(iROI).per = [0, ICA_n(run).area_var(iROI).var'];
           ICA_n(run).area_per(iROI).per = diff(ICA_n(run).area_per(iROI).per)*100;
    
       else
           
           area_ic = zeros(nTR,1);
           
           ICA_n(run).area_var(iROI).var = zeros(nTR,1);
           ICA_n(run).area_per(iROI).per = zeros(nTR,1);
           
       end
       
       ICA_n(run).nICROI(iROI) = size(area_ic,2);
       
       high_area_ic = area_ic(1:nTR,:);
       rest_area_ic = area_ic(nTR+1:end,:);
       
       area_occipital_high = [high_occipital_ic,high_area_ic];
       area_occipital_rest = [rest_occipital_ic,rest_area_ic];
         
       area_parietal_high = [high_parietal_ic,high_area_ic];
       area_parietal_rest = [rest_parietal_ic,rest_area_ic];
        
       area_temporal_high = [high_temporal_ic,high_area_ic];
       area_temporal_rest = [rest_temporal_ic,rest_area_ic];
        
       area_frontal_high = [high_frontal_ic,high_area_ic];
       area_frontal_rest = [rest_frontal_ic,rest_area_ic];
       
       [ICA_FC(run).ROI(iROI).rho_high_occipital,ICA_FC(run).ROI(iROI).pval_high_occipital] = corr(area_occipital_high);
       [ICA_FC(run).ROI(iROI).rho_rest_occipital,ICA_FC(run).ROI(iROI).pval_rest_occipital] = corr(area_occipital_rest);

       [ICA_FC(run).ROI(iROI).rho_high_parietal,ICA_FC(run).ROI(iROI).pval_high_parietal] = corr(area_parietal_high);
       [ICA_FC(run).ROI(iROI).rho_rest_parietal,ICA_FC(run).ROI(iROI).pval_rest_parietal] = corr(area_parietal_rest);

       [ICA_FC(run).ROI(iROI).rho_high_temporal,ICA_FC(run).ROI(iROI).pval_high_temporal] = corr(area_temporal_high);
       [ICA_FC(run).ROI(iROI).rho_rest_temporal,ICA_FC(run).ROI(iROI).pval_rest_temporal] = corr(area_temporal_rest);

       [ICA_FC(run).ROI(iROI).rho_high_frontal,ICA_FC(run).ROI(iROI).pval_high_frontal] = corr(area_frontal_high);
       [ICA_FC(run).ROI(iROI).rho_rest_frontal,ICA_FC(run).ROI(iROI).pval_rest_frontal] = corr(area_frontal_rest);

       clear area_ica
       clear high_area_ic
       clear rest_area_ic
       
    end
  
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-AAL-',preprocessed_step,'.mat'),'ICA_n','ICA_FC','ICA_preprocessed');

end

function plotResults(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for run=1:length(ICA_folder)
    
   max_rho(run) = max(max(ICA_FC(run).rho));
   min_rho(run) = min(min(ICA_FC(run).rho));
    
end

max_all = max(max_rho);
min_all = min(min_rho);

for run=1:length(ICA_folder)
    
   f = figure;
   
   imagesc(ICA_FC(run).rho,[min_all max_all]);
   
   xlabel('occipital - parietal - temporal - frontal');
   ylabel('frontal - temporal - parietal - occipital');
   
   title(ICA_label{run});
   
   colorbar;
   
   hold on;
   
   p1(1) = ICA_n(run).nOccipital;
   p2(1) = 1;
   p1(2) = ICA_n(run).nOccipital;
   p2(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(run).nOccipital;
   p1(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   p2(2) = ICA_n(run).nOccipital;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = ICA_n(run).nOccipital + ICA_n(run).nParietal;
   p2(1) = 1;
   p1(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal;
   p2(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(run).nOccipital + ICA_n(run).nParietal;
   p1(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   p2(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal;
   p2(1) = 1;
   p1(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal;
   p2(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal;
   p1(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal + ICA_n(run).nFrontal;
   p2(2) = ICA_n(run).nOccipital + ICA_n(run).nParietal + ICA_n(run).nTemporal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.pdf'));
   
end

end

function plotCumulative(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-AAL-',preprocessed_step,'.mat'));

nROI = 90;
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for run=1:2
    
    NoA = ICA_n(run).nICOccipital;
    NpA = ICA_n(run).nICParietal;
    NtA = ICA_n(run).nICTemporal;
    NfA = ICA_n(run).nICFrontal;
    
    for iROI=1:nROI
        
        area_label = AAL_ROI(iROI).Nom_L;
        area_label = strrep(area_label,'_','-');
        
        NAA = ICA_n(run).nICROI(iROI);
        
        NA = min( [NoA NpA NtA NfA NAA] ); 
        
%         SOOA = zeros(1,NA);
%         SPPA = zeros(1,NA);
%         STTA = zeros(1,NA);
%         SFFA = zeros(1,NA);

        AOA_high = ICA_FC(run).ROI(iROI).rho_high_occipital(1:NoA,NoA+1:end);
        APA_high = ICA_FC(run).ROI(iROI).rho_high_parietal(1:NpA,NpA+1:end);
        ATA_high = ICA_FC(run).ROI(iROI).rho_high_temporal(1:NtA,NtA+1:end);
        AFA_high = ICA_FC(run).ROI(iROI).rho_high_frontal(1:NfA,NfA+1:end);
        
        AOA_rest = ICA_FC(run).ROI(iROI).rho_rest_occipital(1:NoA,NoA+1:end);
        APA_rest = ICA_FC(run).ROI(iROI).rho_rest_parietal(1:NpA,NpA+1:end);
        ATA_rest = ICA_FC(run).ROI(iROI).rho_rest_temporal(1:NtA,NtA+1:end);
        AFA_rest = ICA_FC(run).ROI(iROI).rho_rest_frontal(1:NfA,NfA+1:end);
        
        SAOA_high = zeros(1,NA);
        SAPA_high = zeros(1,NA);
        SATA_high = zeros(1,NA);
        SAFA_high = zeros(1,NA);
        
        SAOA_rest = zeros(1,NA);
        SAPA_rest = zeros(1,NA);
        SATA_rest = zeros(1,NA);
        SAFA_rest = zeros(1,NA);
        
        for nsc = 2:NA    % number of components to be summed
            
                for ic = 2 : nsc   % loop over components in i-direction
                    
                    for jc = 1 : ic-1   % loop over componnets in j-direction
                        
                        SAOA_high(nsc) = SAOA_high(nsc) + abs( AOA_high(ic,jc) );
                        SAPA_high(nsc) = SAPA_high(nsc) + abs( APA_high(ic,jc) );
                        SATA_high(nsc) = SATA_high(nsc) + abs( ATA_high(ic,jc) );
                        SAFA_high(nsc) = SAFA_high(nsc) + abs( AFA_high(ic,jc) );
                        
                        SAOA_rest(nsc) = SAOA_rest(nsc) + abs( AOA_rest(ic,jc) );
                        SAPA_rest(nsc) = SAPA_rest(nsc) + abs( APA_rest(ic,jc) );
                        SATA_rest(nsc) = SATA_rest(nsc) + abs( ATA_rest(ic,jc) );
                        SAFA_rest(nsc) = SAFA_rest(nsc) + abs( AFA_rest(ic,jc) );
                        
                    end
                    
                end
    
              SAOA_high(nsc) = 2 * SAOA_high(nsc) / (nsc*nsc);   % normalization
              SAPA_high(nsc) = 2 * SAPA_high(nsc) / (nsc*nsc);   % normalization
              SATA_high(nsc) = 2 * SATA_high(nsc) / (nsc*nsc);   % normalization
              SAFA_high(nsc) = 2 * SAFA_high(nsc) / (nsc*nsc);   % normalization
              
              SAOA_rest(nsc) = 2 * SAOA_rest(nsc) / (nsc*nsc);   % normalization
              SAPA_rest(nsc) = 2 * SAPA_rest(nsc) / (nsc*nsc);   % normalization
              SATA_rest(nsc) = 2 * SATA_rest(nsc) / (nsc*nsc);   % normalization
              SAFA_rest(nsc) = 2 * SAFA_rest(nsc) / (nsc*nsc);   % normalization

%             sao_high = abs( AOA_high(1:nsc,1:nsc) );
%             sap_high = abs( APA_high(1:nsc,1:nsc) );
%             sat_high = abs( ATA_high(1:nsc,1:nsc) );
%             saf_high = abs( AFA_high(1:nsc,1:nsc) );
%             
%             sao_rest = abs( AOA_rest(1:nsc,1:nsc) );
%             sap_rest = abs( APA_rest(1:nsc,1:nsc) );
%             sat_rest = abs( ATA_rest(1:nsc,1:nsc) );
%             saf_rest = abs( AFA_rest(1:nsc,1:nsc) );
%             
%             SAOA_high(nsc) = ( sum(sao_high(:)) ) / (nsc*nsc);
%             SAPA_high(nsc) = ( sum(sap_high(:)) ) / (nsc*nsc);
%             SATA_high(nsc) = ( sum(sat_high(:)) ) / (nsc*nsc);
%             SAFA_high(nsc) = ( sum(saf_high(:)) ) / (nsc*nsc);
% 
%             SAOA_rest(nsc) = ( sum(sao_rest(:)) ) / (nsc*nsc);
%             SAPA_rest(nsc) = ( sum(sap_rest(:)) ) / (nsc*nsc);
%             SATA_rest(nsc) = ( sum(sat_rest(:)) ) / (nsc*nsc);
%             SAFA_rest(nsc) = ( sum(saf_rest(:)) ) / (nsc*nsc);

        end
    
    f = figure;
    
    fs = 6;
    
    subplot(1,4,1);
    
    hold on;
    plot(1:NA, SAOA_high, 'ko', 'LineWidth', 2 );
    plot(1:NA, SAOA_rest, 'bo', 'LineWidth', 2 );
    hold off;
    
    h = legend('High Attention', 'Resting State');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(area_label,'-','Occipital Lobe'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    subplot(1,4,2);
    
    hold on;
    plot(1:NA, SAPA_high, 'ko', 'LineWidth', 2 );
    plot(1:NA, SAPA_rest, 'bo', 'LineWidth', 2 );
    hold off;
    
    h = legend('High Attention', 'Resting State');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(area_label,'-','Parietal Lobe'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );
    
    subplot(1,4,3);
    
    hold on;
    plot(1:NA, SATA_high, 'ko', 'LineWidth', 2 );
    plot(1:NA, SATA_rest, 'bo', 'LineWidth', 2 );
    hold off;
    
    h = legend('High Attention', 'Resting State');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(area_label,'-','Temporal Lobe'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );
    
    subplot(1,4,4);
    
    hold on;
    plot(1:NA, SAFA_high, 'ko', 'LineWidth', 2 );
    plot(1:NA, SAFA_rest, 'bo', 'LineWidth', 2 );
    hold off;
    
    h = legend('High Attention', 'Resting State');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(area_label,'-','Frontal Lobe'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes-AAL','-',area_label,'-run-',int2str(run),'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes-AAL','-',area_label,'-run-',int2str(run),'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes-AAL','-',area_label,'-run-',int2str(run),'-',preprocessed_step,'.pdf'));
   
    close all;
    
    end
    
end
   
end

function plotSixCumulativeTwoConditions(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for run=1:2
    
    NoA = ICA_n(run).nOccipital;
    NpA = ICA_n(run).nParietal;
    NtA = ICA_n(run).nTemporal;
    NfA = ICA_n(run).nFrontal;
    
    NA_all(run) = min( [NoA NpA NtA NfA] ); 
    
end

%NA = min(NA_all);

for run=1:2
    
    NoA = ICA_n(run).nOccipital;
    NpA = ICA_n(run).nParietal;
    NtA = ICA_n(run).nTemporal;
    NfA = ICA_n(run).nFrontal;
    
    NA = NA_all(run);
   
    %%% HIGH 
    
    Run(run).Condition(1).OOA = ICA_FC(run).rho_high(1:NoA,1:NoA);
    Run(run).Condition(1).OPA = ICA_FC(run).rho_high(1:NoA,NoA+1:NoA+NpA);
    Run(run).Condition(1).OTA = ICA_FC(run).rho_high(1:NoA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(1).OFA = ICA_FC(run).rho_high(1:NoA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(1).PPA = ICA_FC(run).rho_high(NoA+1:NoA+NpA,NoA+1:NoA+NpA);
    Run(run).Condition(1).PTA = ICA_FC(run).rho_high(NoA+1:NoA+NpA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(1).PFA = ICA_FC(run).rho_high(NoA+1:NoA+NpA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(1).TTA = ICA_FC(run).rho_high(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(1).TFA = ICA_FC(run).rho_high(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(1).FFA = ICA_FC(run).rho_high(NoA+NpA+NtA+1:NoA+NpA+NtA+NfA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA); 

    %%% REST
    
    Run(run).Condition(2).OOA = ICA_FC(run).rho_rest(1:NoA,1:NoA);
    Run(run).Condition(2).OPA = ICA_FC(run).rho_rest(1:NoA,NoA+1:NoA+NpA);
    Run(run).Condition(2).OTA = ICA_FC(run).rho_rest(1:NoA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(2).OFA = ICA_FC(run).rho_rest(1:NoA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(2).PPA = ICA_FC(run).rho_rest(NoA+1:NoA+NpA,NoA+1:NoA+NpA);
    Run(run).Condition(2).PTA = ICA_FC(run).rho_rest(NoA+1:NoA+NpA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(2).PFA = ICA_FC(run).rho_rest(NoA+1:NoA+NpA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(2).TTA = ICA_FC(run).rho_rest(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+1:NoA+NpA+NtA);
    Run(run).Condition(2).TFA = ICA_FC(run).rho_rest(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Run(run).Condition(2).FFA = ICA_FC(run).rho_rest(NoA+NpA+NtA+1:NoA+NpA+NtA+NfA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA); 

    for iCondition=1:2
        
        Run(run).Condition(iCondition).SOOA = zeros(1,NA);
        Run(run).Condition(iCondition).SOPA = zeros(1,NA);
        Run(run).Condition(iCondition).SOTA = zeros(1,NA);
        Run(run).Condition(iCondition).SOFA = zeros(1,NA);
        Run(run).Condition(iCondition).SPPA = zeros(1,NA);
        Run(run).Condition(iCondition).SPTA = zeros(1,NA);
        Run(run).Condition(iCondition).SPFA = zeros(1,NA);
        Run(run).Condition(iCondition).STTA = zeros(1,NA);
        Run(run).Condition(iCondition).STFA = zeros(1,NA);
        Run(run).Condition(iCondition).SFFA = zeros(1,NA);

        for nsc = 2:NA    % number of components to be summed

%             for ic = 2 : nsc   % loop over components in i-direction
%                 for jc = 1 : ic-1   % loop over componnets in j-direction
% 
%                     Run(run).Condition(iCondition).SOOA(nsc) = Run(run).Condition(iCondition).SOOA(nsc) + abs(  Run(run).Condition(iCondition).OOA(ic,jc) );
%                     Run(run).Condition(iCondition).SOPA(nsc) = Run(run).Condition(iCondition).SOPA(nsc) + abs(  Run(run).Condition(iCondition).OPA(ic,jc) );
%                     Run(run).Condition(iCondition).SOTA(nsc) = Run(run).Condition(iCondition).SOTA(nsc) + abs(  Run(run).Condition(iCondition).OTA(ic,jc) );
%                     Run(run).Condition(iCondition).SOFA(nsc) = Run(run).Condition(iCondition).SOFA(nsc) + abs(  Run(run).Condition(iCondition).OFA(ic,jc) );
%                     Run(run).Condition(iCondition).SPPA(nsc) = Run(run).Condition(iCondition).SPPA(nsc) + abs(  Run(run).Condition(iCondition).PPA(ic,jc) );
%                     Run(run).Condition(iCondition).SPTA(nsc) = Run(run).Condition(iCondition).SPTA(nsc) + abs(  Run(run).Condition(iCondition).PTA(ic,jc) );
%                     Run(run).Condition(iCondition).SPFA(nsc) = Run(run).Condition(iCondition).SPFA(nsc) + abs(  Run(run).Condition(iCondition).PFA(ic,jc) );
%                     Run(run).Condition(iCondition).STTA(nsc) = Run(run).Condition(iCondition).STTA(nsc) + abs(  Run(run).Condition(iCondition).TTA(ic,jc) );
%                     Run(run).Condition(iCondition).STFA(nsc) = Run(run).Condition(iCondition).STFA(nsc) + abs(  Run(run).Condition(iCondition).TFA(ic,jc) );
%                     Run(run).Condition(iCondition).SFFA(nsc) = Run(run).Condition(iCondition).SFFA(nsc) + abs(  Run(run).Condition(iCondition).FFA(ic,jc) );
% 
%                 end
%             end

            soo = abs(  Run(run).Condition(iCondition).OOA(1:nsc,1:nsc) );
            sop = abs(  Run(run).Condition(iCondition).OPA(1:nsc,1:nsc) );
            sot = abs(  Run(run).Condition(iCondition).OTA(1:nsc,1:nsc) );
            sof = abs(  Run(run).Condition(iCondition).OFA(1:nsc,1:nsc) );
            spp = abs(  Run(run).Condition(iCondition).PPA(1:nsc,1:nsc) );
            spt = abs(  Run(run).Condition(iCondition).PTA(1:nsc,1:nsc) );
            spf = abs(  Run(run).Condition(iCondition).PFA(1:nsc,1:nsc) );
            stt = abs(  Run(run).Condition(iCondition).TTA(1:nsc,1:nsc) );
            stf = abs(  Run(run).Condition(iCondition).TFA(1:nsc,1:nsc) );
            sff = abs(  Run(run).Condition(iCondition).FFA(1:nsc,1:nsc) );

            Run(run).Condition(iCondition).SOOA(nsc) = ( sum(soo(:)) - sum(diag(soo)) ) / (nsc*nsc-nsc);
            Run(run).Condition(iCondition).SOPA(nsc) = ( sum(sop(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).SOTA(nsc) = ( sum(sot(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).SOFA(nsc) = ( sum(sof(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).SPPA(nsc) = ( sum(spp(:)) - sum(diag(spp)) ) / (nsc*nsc-nsc);
            Run(run).Condition(iCondition).SPTA(nsc) = ( sum(spt(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).SPFA(nsc) = ( sum(spf(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).STTA(nsc) = ( sum(stt(:)) - sum(diag(stt)) ) / (nsc*nsc-nsc);
            Run(run).Condition(iCondition).STFA(nsc) = ( sum(stf(:))                  ) / (nsc*nsc);
            Run(run).Condition(iCondition).SFFA(nsc) = ( sum(sff(:)) - sum(diag(sff)) ) / (nsc*nsc-nsc);
            
%             Run(run).Condition(iCondition).SOOA(nsc) = 2 * Run(run).Condition(iCondition).SOOA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SOPA(nsc) = 2 * Run(run).Condition(iCondition).SOPA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SOTA(nsc) = 2 * Run(run).Condition(iCondition).SOTA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SOFA(nsc) = 2 * Run(run).Condition(iCondition).SOFA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SPPA(nsc) = 2 * Run(run).Condition(iCondition).SPPA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SPTA(nsc) = 2 * Run(run).Condition(iCondition).SPTA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SPFA(nsc) = 2 * Run(run).Condition(iCondition).SPFA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).STTA(nsc) = 2 * Run(run).Condition(iCondition).STTA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).STFA(nsc) = 2 * Run(run).Condition(iCondition).STFA(nsc) / (nsc*nsc);   % normalization
%             Run(run).Condition(iCondition).SFFA(nsc) = 2 * Run(run).Condition(iCondition).SFFA(nsc) / (nsc*nsc);   % normalization
       
        end
    
    end

end
    
run = 1;
plotTwoCumulative(Run(run).Condition(1).SOOA,Run(run).Condition(2).SOOA,'HighAttention','RestingState','O-O',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SPPA,Run(run).Condition(2).SPPA,'HighAttention','RestingState','P-P',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SFFA,Run(run).Condition(2).SFFA,'HighAttention','RestingState','F-F',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SOPA,Run(run).Condition(2).SOPA,'HighAttention','RestingState','O-P',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SPFA,Run(run).Condition(2).SPFA,'HighAttention','RestingState','P-F',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SOFA,Run(run).Condition(2).SOFA,'HighAttention','RestingState','O-F',NA_all(run),run,settings,preprocessed_step);

run = 2;
plotTwoCumulative(Run(run).Condition(1).SOOA,Run(run).Condition(2).SOOA,'HighAttention','RestingState','O-O',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SPPA,Run(run).Condition(2).SPPA,'HighAttention','RestingState','P-P',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SFFA,Run(run).Condition(2).SFFA,'HighAttention','RestingState','F-F',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SOPA,Run(run).Condition(2).SOPA,'HighAttention','RestingState','O-P',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SPFA,Run(run).Condition(2).SPFA,'HighAttention','RestingState','P-F',NA_all(run),run,settings,preprocessed_step);

plotTwoCumulative(Run(run).Condition(1).SOFA,Run(run).Condition(2).SOFA,'HighAttention','RestingState','O-F',NA_all(run),run,settings,preprocessed_step);

  
end

function plotTwoCumulative(Cumulative1,Cumulative2,label1,label2,titlelabel,NA,run,settings,preprocessed_step)

    fs = 14;
    
    f = figure;
    
    hold on;
    plot(1:NA, Cumulative1, 'ko', 'LineWidth', 2 );
    plot(1:NA, Cumulative2, 'bo', 'LineWidth', 2 );
    
    h = legend(label1, label2);
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(titlelabel, 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-groupICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-groupICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-groupICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.pdf'));
   

end

function plotVarianceHistogram(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-AAL-',preprocessed_step,'.mat'));

load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

nROI = 90;

for iRun=1:2

    %%% FOR ROIs
    
    %for iROI=1:nROI
        for iROI=65:68
        
        area_label = AAL_ROI(iROI).Nom_L;
        area_label = strrep(area_label,'_','-');
    
        nComponents = ICA_n(iRun).nICROI(iROI);
        
        IC_percent_variance = ICA_n(iRun).area_per(iROI).per(1,1:nComponents);
        
        M = 1000;
        
        f = figure;
        
        hist(IC_percent_variance,M);
        
        xlabel('Percentage of Variance Explained');
        ylabel('Number of Components');
        title(strcat(area_label,'-','run:',int2str(iRun)));
        %xlim([0 1.5]);
        
        print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','Hist-Percent-Variance-',preprocessed_step,'.jpeg'));
        print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','Hist-Percent-Variance-',preprocessed_step,'.eps'));
        
        close all;
        
    end
    
    %%% FOR LOBES
    
    %%% Occipital
    nComponents = ICA_n(iRun).nICOccipital;
    IC_percent_variance = ICA_n(iRun).occipital_per(1,1:nComponents);
    
    M = 1000;
        
    f = figure;

    hist(IC_percent_variance,M);

    xlabel('Percentage of Variance Explained');
    ylabel('Number of Components');
    title(strcat('Occipital','-','run:',int2str(iRun)));
    %xlim([0 1.5]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Occipital','-','Hist-Percent-Variance-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Occipital','-','Hist-Percent-Variance-',preprocessed_step,'.eps'));

    %%% Parietal
    nComponents = ICA_n(iRun).nICParietal;
    IC_percent_variance = ICA_n(iRun).parietal_per(1,1:nComponents);
    
    M = 1000;
        
    f = figure;

    hist(IC_percent_variance,M);

    xlabel('Percentage of Variance Explained');
    ylabel('Number of Components');
    title(strcat('Parietal','-','run:',int2str(iRun)));
    %xlim([0 1.5]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Parietal','-','Hist-Percent-Variance-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Parietal','-','Hist-Percent-Variance-',preprocessed_step,'.eps'));

    %%% Temporal
    nComponents = ICA_n(iRun).nICTemporal;
    IC_percent_variance = ICA_n(iRun).temporal_per(1,1:nComponents);
    
    M = 1000;
        
    f = figure;

    hist(IC_percent_variance,M);

    xlabel('Percentage of Variance Explained');
    ylabel('Number of Components');
    title(strcat('Temporal','-','run:',int2str(iRun)));
    %xlim([0 1.5]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Temporal','-','Hist-Percent-Variance-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Temporal','-','Hist-Percent-Variance-',preprocessed_step,'.eps'));

    %%% Frontal
    nComponents = ICA_n(iRun).nICFrontal;
    IC_percent_variance = ICA_n(iRun).frontal_per(1,1:nComponents);
    
    M = 1000;
        
    f = figure;

    hist(IC_percent_variance,M);

    xlabel('Percentage of Variance Explained');
    ylabel('Number of Components');
    title(strcat('Frontal','-','run:',int2str(iRun)));
    %xlim([0 1.5]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Frontal','-','Hist-Percent-Variance-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-','Frontal','-','Hist-Percent-Variance-',preprocessed_step,'.eps'));

    close all;
    
end


end

function plotTimeCourses(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-AAL-',preprocessed_step,'.mat'));

load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ICA_preprocessed = 'groupICA';

nTR = 300;

area_lobe_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area_lobe');
area_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area');

nROI = 90;

nIC = 40;

for iRun=1:2

    %%% FOR ROIs
    
    %for iROI=1:nROI
        for iROI=65:68
            
            area_label = AAL_ROI(iROI).Nom_L;
            area_label = strrep(area_label,'_','-');
            
            nICOccipital = ICA_n(iRun).nICOccipital;
            
            for iIC=1:nICOccipital
                
                sum_rho_high(iIC) = sum(ICA_FC(iRun).ROI(iROI).rho_high_occipital(iIC,nICOccipital+1:end));
                sum_rho_rest(iIC) = sum(ICA_FC(iRun).ROI(iROI).rho_rest_occipital(iIC,nICOccipital+1:end));
            
            end
            
            [Y, I] = sort(sum_rho_high);
            max_rho_high = I(end);
            
            clear Y;
            clear I;
            
            [Y, I] = sort(sum_rho_rest);
            max_rho_rest = I(end);
            
            clear Y;
            clear I;
            
            clear sum_rho_high;
            clear sum_rho_rest;
            
            occipital_ic = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(iRun),'-','occipital_lobe','\','melodic_mix'));
            
            time_course_high = occipital_ic(1:nTR,max_rho_high);
            time_course_rest = occipital_ic(nTR+1:end,max_rho_rest);
            
            [Y, I] = sort(ICA_FC(iRun).ROI(iROI).rho_high_occipital(max_rho_high,nICOccipital+1:end));
            area_ic_high_idx = I(end-9:end);
            
            clear Y;
            clear I;
            
            [Y, I] = sort(ICA_FC(iRun).ROI(iROI).rho_rest_occipital(max_rho_rest,nICOccipital+1:end));
            area_ic_rest_idx = I(end-9:end);
            
            clear Y;
            clear I;
            
            area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(iRun),'-','aal2std_',int2str(iROI),'\','melodic_mix'));
             
            high_area_ic = area_ic(1:nTR,:);
            rest_area_ic = area_ic(nTR+1:end,:);
            
            f = figure;
            
            fs = 6;
            
            for iIC=1:length(area_ic_high_idx)
               
                subplot(5,2,iIC);
                
                plot(time_course_high,'r');
                hold on;
                plot(high_area_ic(:,area_ic_high_idx(iIC)),'b');
                
                title(strcat(area_label,'-IC:',int2str(area_ic_high_idx(iIC)),'-Occipital','-IC:',int2str(max_rho_high),'-High','-run:',int2str(iRun)),'FontSize',fs);
                xlabel('TRs','FontSize',fs);
                ylabel('BOLD','FontSize',fs);
                
            end
            
            print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-High-',preprocessed_step,'.jpeg'));
            print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-High-',preprocessed_step,'.eps'));
            print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-High-',preprocessed_step,'.pdf'));
            
            close all;
            
            f = figure;
            
            fs = 6;
            
            for iIC=1:length(area_ic_rest_idx)
               
                subplot(5,2,iIC);
                
                plot(time_course_rest,'r');
                hold on;
                plot(high_area_ic(:,area_ic_rest_idx(iIC)),'b');
                
                title(strcat(area_label,'-IC:',int2str(area_ic_rest_idx(iIC)),'-Occipital','-IC:',int2str(max_rho_rest),'-Rest','-run:',int2str(iRun)),'FontSize',fs);
                xlabel('TRs','FontSize',fs);
                ylabel('BOLD','FontSize',fs);
                
            end
            
            print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-Rest-',preprocessed_step,'.jpeg'));
            print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-Rest-',preprocessed_step,'.eps'));
            print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-groupICA-Four-lobes-AAL','-','Run-',int2str(iRun),'-',area_label,'-','IC-Time-Course-Occipital-Rest-',preprocessed_step,'.pdf'));
            
            close all;
            
        end
        
end

end