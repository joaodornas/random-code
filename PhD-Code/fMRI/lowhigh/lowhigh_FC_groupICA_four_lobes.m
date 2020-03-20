function lowhigh_FC_groupICA_four_lobes


%settings_jan_0805;
settings_elena_2905;

preprocessed_step = 'residual';
%preprocessed_step = 'filtered';

%doTheMath(settings,preprocessed_step);

%plotResults(settings,preprocessed_step);

%plotCumulative(settings,preprocessed_step);

plotSixCumulativeTwoConditions(settings,preprocessed_step);

end

function doTheMath(settings,preprocessed_step)

ICA_preprocessed = 'groupICA';

nTR = 300;

preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area_lobe');

for run=1:2
    
    disp(strcat('run:',int2str(run)));
    
    occipital_ic = load(strcat(preprocessed_folder,'\','high-rest-run-',int2str(run),'-','occipital_lobe','\','melodic_mix'));
    parietal_ic = load(strcat(preprocessed_folder,'\','high-rest-run-',int2str(run),'-','parietal_lobe','\','melodic_mix'));
    temporal_ic = load(strcat(preprocessed_folder,'\','high-rest-run-',int2str(run),'-','temporal_lobe','\','melodic_mix'));
    frontal_ic = load(strcat(preprocessed_folder,'\','high-rest-run-',int2str(run),'-','frontal_lobe','\','melodic_mix'));
    
    ICA_n(run).nOccipital = size(occipital_ic,2);
    ICA_n(run).nParietal = size(parietal_ic,2);
    ICA_n(run).nTemporal = size(temporal_ic,2);
    ICA_n(run).nFrontal = size(frontal_ic,2);
    
    high_occipital_ic = occipital_ic(1:nTR,:);
    high_parietal_ic = parietal_ic(1:nTR,:);
    high_temporal_ic = temporal_ic(1:nTR,:);
    high_frontal_ic = frontal_ic(1:nTR,:);
    
    rest_occipital_ic = occipital_ic(nTR+1:end,:);
    rest_parietal_ic = parietal_ic(nTR+1:end,:);
    rest_temporal_ic = temporal_ic(nTR+1:end,:);
    rest_frontal_ic = frontal_ic(nTR+1:end,:);
    
    allComponents_high = [high_occipital_ic,high_parietal_ic,high_temporal_ic,high_frontal_ic];
    allComponents_rest = [rest_occipital_ic,rest_parietal_ic,rest_temporal_ic,rest_frontal_ic];
    
    [ICA_FC(run).rho_high,ICA_FC(run).pval_high] = corr(allComponents_high);
    [ICA_FC(run).rho_rest,ICA_FC(run).pval_rest] = corr(allComponents_rest);
  
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'),'ICA_n','ICA_FC','ICA_preprocessed');

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

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for run=1:6
    
    NoA = ICA_n(run).nOccipital;
    NpA = ICA_n(run).nParietal;
    NtA = ICA_n(run).nTemporal;
    NfA = ICA_n(run).nFrontal;
   
    OOA = ICA_FC(run).rho(1:NoA,1:NoA);
    OPA = ICA_FC(run).rho(1:NoA,NoA+1:NoA+NpA);
    OTA = ICA_FC(run).rho(1:NoA,NoA+NpA+1:NoA+NpA+NtA);
    OFA = ICA_FC(run).rho(1:NoA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    PPA = ICA_FC(run).rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);
    PTA = ICA_FC(run).rho(NoA+1:NoA+NpA,NoA+NpA+1:NoA+NpA+NtA);
    PFA = ICA_FC(run).rho(NoA+1:NoA+NpA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    TTA = ICA_FC(run).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+1:NoA+NpA+NtA);
    TFA = ICA_FC(run).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    FFA = ICA_FC(run).rho(NoA+NpA+NtA+1:NoA+NpA+NtA+NfA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA); 
    
    NA = min( [NoA NpA NtA NfA] );  

    SOOA = zeros(1,NA);
    SOPA = zeros(1,NA);
    SOTA = zeros(1,NA);
    SOFA = zeros(1,NA);
    SPPA = zeros(1,NA);
    SPTA = zeros(1,NA);
    SPFA = zeros(1,NA);
    STTA = zeros(1,NA);
    STFA = zeros(1,NA);
    SFFA = zeros(1,NA);

    for nsc = 2:NA    % number of components to be summed

        for ic = 2 : nsc   % loop over components in i-direction
            for jc = 1 : ic-1   % loop over componnets in j-direction
                
                SOOA(nsc) = SOOA(nsc) + abs( OOA(ic,jc) );
                SOPA(nsc) = SOPA(nsc) + abs( OPA(ic,jc) );
                SOTA(nsc) = SOTA(nsc) + abs( OTA(ic,jc) );
                SOFA(nsc) = SOFA(nsc) + abs( OFA(ic,jc) );
                SPPA(nsc) = SPPA(nsc) + abs( PPA(ic,jc) );
                SPTA(nsc) = SPTA(nsc) + abs( PTA(ic,jc) );
                SPFA(nsc) = SPFA(nsc) + abs( PFA(ic,jc) );
                STTA(nsc) = STTA(nsc) + abs( TTA(ic,jc) );
                STFA(nsc) = STFA(nsc) + abs( TFA(ic,jc) );
                SFFA(nsc) = SFFA(nsc) + abs( FFA(ic,jc) );

            end
        end
        
        SOOA(nsc) = 2 * SOOA(nsc) / (nsc*nsc);   % normalization
        SOPA(nsc) = 2 * SOPA(nsc) / (nsc*nsc);   % normalization
        SOTA(nsc) = 2 * SOTA(nsc) / (nsc*nsc);   % normalization
        SOFA(nsc) = 2 * SOFA(nsc) / (nsc*nsc);   % normalization
        SPPA(nsc) = 2 * SPPA(nsc) / (nsc*nsc);   % normalization
        SPTA(nsc) = 2 * SPTA(nsc) / (nsc*nsc);   % normalization
        SPFA(nsc) = 2 * SPFA(nsc) / (nsc*nsc);   % normalization
        STTA(nsc) = 2 * STTA(nsc) / (nsc*nsc);   % normalization
        STFA(nsc) = 2 * STFA(nsc) / (nsc*nsc);   % normalization
        SFFA(nsc) = 2 * SFFA(nsc) / (nsc*nsc);   % normalization
end
    
    f = figure;
    
    fs = 14;
    
    subplot(2,2,1);
    
    hold on;
    plot(1:NA, SOOA, 'ko', 'LineWidth', 2 );
    plot(1:NA, SOPA, 'bo', 'LineWidth', 2 );
    plot(1:NA, SOTA, 'ro', 'LineWidth', 2 );
    plot(1:NA, SOFA, 'yo', 'LineWidth', 2 );
    hold off;
    
    h = legend('O-O', 'O-P', 'O-T', 'O-F');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(ICA_label{run},'-','Occipital'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    subplot(2,2,2);
    
    hold on;
    plot(1:NA, SPPA, 'ko', 'LineWidth', 2 );
    plot(1:NA, SOPA, 'bo', 'LineWidth', 2 );
    plot(1:NA, SPTA, 'ro', 'LineWidth', 2 );
    plot(1:NA, SPFA, 'yo', 'LineWidth', 2 );
    hold off;
    
    h = legend('P-P', 'P-O', 'P-T', 'P-F');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(ICA_label{run},'-','Parietal'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    subplot(2,2,3);
    
    hold on;
    plot(1:NA, STTA, 'ko', 'LineWidth', 2 );
    plot(1:NA, SOTA, 'bo', 'LineWidth', 2 );
    plot(1:NA, SPTA, 'ro', 'LineWidth', 2 );
    plot(1:NA, STFA, 'yo', 'LineWidth', 2 );
    hold off;
    
    h = legend('T-T', 'T-O', 'T-P', 'T-F');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(ICA_label{run},'-','Temporal'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    subplot(2,2,4);
    
    hold on;
    plot(1:NA, SFFA, 'ko', 'LineWidth', 2 );
    plot(1:NA, SOFA, 'bo', 'LineWidth', 2 );
    plot(1:NA, SPFA, 'ro', 'LineWidth', 2 );
    plot(1:NA, STFA, 'yo', 'LineWidth', 2 );
    hold off;
    
    h = legend('F-F', 'F-O', 'F-P', 'F-T');
    
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(strcat(ICA_label{run},'-','Frontal'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{run},'-',preprocessed_step,'.pdf'));
   
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
