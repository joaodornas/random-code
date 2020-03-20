function lowhigh_FC_ICA_four_lobes


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

preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed);

ICA_preprocessed = 'ICA_lobe';
%ICA_preprocessed = 'ICA_lobe_filtered';

ICA_folder{1} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(1).folder,'\','FSL','\',ICA_preprocessed);
ICA_folder{2} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(2).folder,'\','FSL','\',ICA_preprocessed);

ICA_folder{3} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(1).folder,'\','FSL','\',ICA_preprocessed);
ICA_folder{4} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(2).folder,'\','FSL','\',ICA_preprocessed);

ICA_folder{5} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(1).folder,'\','FSL','\',ICA_preprocessed);
ICA_folder{6} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(2).folder,'\','FSL','\',ICA_preprocessed);

ICA_label{1} = 'High-Attention-Run-1';
ICA_label{2} = 'High-Attention-Run-2';
ICA_label{3} = 'Low-Attention-Run-1';
ICA_label{4} = 'Low-Attention-Run-2';
ICA_label{5} = 'RestingState-Run-1';
ICA_label{6} = 'RestingState-Run-2';

for ifolder=1:length(ICA_folder)
    
    disp(ICA_label{ifolder});
    
    occipital_ic = load(strcat(ICA_folder{ifolder},'\','occipital_lobe','\','melodic_mix'));
    parietal_ic = load(strcat(ICA_folder{ifolder},'\','parietal_lobe','\','melodic_mix'));
    temporal_ic = load(strcat(ICA_folder{ifolder},'\','temporal_lobe','\','melodic_mix'));
    frontal_ic = load(strcat(ICA_folder{ifolder},'\','frontal_lobe','\','melodic_mix'));
    
    ICA_n(ifolder).nOccipital = size(occipital_ic,2);
    ICA_n(ifolder).nParietal = size(parietal_ic,2);
    ICA_n(ifolder).nTemporal = size(temporal_ic,2);
    ICA_n(ifolder).nFrontal = size(frontal_ic,2);
    
    allComponents = [occipital_ic,parietal_ic,temporal_ic,frontal_ic];
    
    [ICA_FC(ifolder).rho,ICA_FC(ifolder).pval] = corr(allComponents);
  
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'),'ICA_label','ICA_folder','ICA_n','ICA_FC','ICA_preprocessed');

end


function plotResults(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for ifolder=1:length(ICA_folder)
    
   max_rho(ifolder) = max(max(ICA_FC(ifolder).rho));
   min_rho(ifolder) = min(min(ICA_FC(ifolder).rho));
    
end

max_all = max(max_rho);
min_all = min(min_rho);

for ifolder=1:length(ICA_folder)
    
   f = figure;
   
   imagesc(ICA_FC(ifolder).rho,[min_all max_all]);
   
   xlabel('occipital - parietal - temporal - frontal');
   ylabel('frontal - temporal - parietal - occipital');
   
   title(ICA_label{ifolder});
   
   colorbar;
   
   hold on;
   
   p1(1) = ICA_n(ifolder).nOccipital;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nOccipital;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nOccipital;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   p2(2) = ICA_n(ifolder).nOccipital;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal + ICA_n(ifolder).nFrontal;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal + ICA_n(ifolder).nTemporal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.pdf'));
   
end

end


function plotCumulative(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for ifolder=1:6
    
    NoA = ICA_n(ifolder).nOccipital;
    NpA = ICA_n(ifolder).nParietal;
    NtA = ICA_n(ifolder).nTemporal;
    NfA = ICA_n(ifolder).nFrontal;
   
    OOA = ICA_FC(ifolder).rho(1:NoA,1:NoA);
    OPA = ICA_FC(ifolder).rho(1:NoA,NoA+1:NoA+NpA);
    OTA = ICA_FC(ifolder).rho(1:NoA,NoA+NpA+1:NoA+NpA+NtA);
    OFA = ICA_FC(ifolder).rho(1:NoA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    PPA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);
    PTA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+NpA+1:NoA+NpA+NtA);
    PFA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    TTA = ICA_FC(ifolder).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+1:NoA+NpA+NtA);
    TFA = ICA_FC(ifolder).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    FFA = ICA_FC(ifolder).rho(NoA+NpA+NtA+1:NoA+NpA+NtA+NfA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA); 
    
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
    title(strcat(ICA_label{ifolder},'-','Occipital'), 'FontSize', fs);
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
    title(strcat(ICA_label{ifolder},'-','Parietal'), 'FontSize', fs);
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
    title(strcat(ICA_label{ifolder},'-','Temporal'), 'FontSize', fs);
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
    title(strcat(ICA_label{ifolder},'-','Frontal'), 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',ICA_label{ifolder},'-',preprocessed_step,'.pdf'));
   
end
   
end

function plotSixCumulativeTwoConditions(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Four-lobes-All-Runs-',preprocessed_step,'.mat'));

for ifolder=1:6
    
    NoA = ICA_n(ifolder).nOccipital;
    NpA = ICA_n(ifolder).nParietal;
    NtA = ICA_n(ifolder).nTemporal;
    NfA = ICA_n(ifolder).nFrontal;
    
    NA_all(ifolder) = min( [NoA NpA NtA NfA] ); 
    
end

NA = min(NA_all);

for ifolder=1:6
    
    NoA = ICA_n(ifolder).nOccipital;
    NpA = ICA_n(ifolder).nParietal;
    NtA = ICA_n(ifolder).nTemporal;
    NfA = ICA_n(ifolder).nFrontal;
   
    Condition(ifolder).OOA = ICA_FC(ifolder).rho(1:NoA,1:NoA);
    Condition(ifolder).OPA = ICA_FC(ifolder).rho(1:NoA,NoA+1:NoA+NpA);
    Condition(ifolder).OTA = ICA_FC(ifolder).rho(1:NoA,NoA+NpA+1:NoA+NpA+NtA);
    Condition(ifolder).OFA = ICA_FC(ifolder).rho(1:NoA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Condition(ifolder).PPA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);
    Condition(ifolder).PTA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+NpA+1:NoA+NpA+NtA);
    Condition(ifolder).PFA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Condition(ifolder).TTA = ICA_FC(ifolder).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+1:NoA+NpA+NtA);
    Condition(ifolder).TFA = ICA_FC(ifolder).rho(NoA+NpA+1:NoA+NpA+NtA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA);
    
    Condition(ifolder).FFA = ICA_FC(ifolder).rho(NoA+NpA+NtA+1:NoA+NpA+NtA+NfA,NoA+NpA+NtA+1:NoA+NpA+NtA+NfA); 

    Condition(ifolder).SOOA = zeros(1,NA);
    Condition(ifolder).SOPA = zeros(1,NA);
    Condition(ifolder).SOTA = zeros(1,NA);
    Condition(ifolder).SOFA = zeros(1,NA);
    Condition(ifolder).SPPA = zeros(1,NA);
    Condition(ifolder).SPTA = zeros(1,NA);
    Condition(ifolder).SPFA = zeros(1,NA);
    Condition(ifolder).STTA = zeros(1,NA);
    Condition(ifolder).STFA = zeros(1,NA);
    Condition(ifolder).SFFA = zeros(1,NA);

    for nsc = 2:NA    % number of components to be summed

        for ic = 2 : nsc   % loop over components in i-direction
            for jc = 1 : ic-1   % loop over componnets in j-direction
                
                Condition(ifolder).SOOA(nsc) = Condition(ifolder).SOOA(nsc) + abs(  Condition(ifolder).OOA(ic,jc) );
                Condition(ifolder).SOPA(nsc) = Condition(ifolder).SOPA(nsc) + abs(  Condition(ifolder).OPA(ic,jc) );
                Condition(ifolder).SOTA(nsc) = Condition(ifolder).SOTA(nsc) + abs(  Condition(ifolder).OTA(ic,jc) );
                Condition(ifolder).SOFA(nsc) = Condition(ifolder).SOFA(nsc) + abs(  Condition(ifolder).OFA(ic,jc) );
                Condition(ifolder).SPPA(nsc) = Condition(ifolder).SPPA(nsc) + abs(  Condition(ifolder).PPA(ic,jc) );
                Condition(ifolder).SPTA(nsc) = Condition(ifolder).SPTA(nsc) + abs(  Condition(ifolder).PTA(ic,jc) );
                Condition(ifolder).SPFA(nsc) = Condition(ifolder).SPFA(nsc) + abs(  Condition(ifolder).PFA(ic,jc) );
                Condition(ifolder).STTA(nsc) = Condition(ifolder).STTA(nsc) + abs(  Condition(ifolder).TTA(ic,jc) );
                Condition(ifolder).STFA(nsc) = Condition(ifolder).STFA(nsc) + abs(  Condition(ifolder).TFA(ic,jc) );
                Condition(ifolder).SFFA(nsc) = Condition(ifolder).SFFA(nsc) + abs(  Condition(ifolder).FFA(ic,jc) );

            end
        end
        
        Condition(ifolder).SOOA(nsc) = 2 * Condition(ifolder).SOOA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SOPA(nsc) = 2 * Condition(ifolder).SOPA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SOTA(nsc) = 2 * Condition(ifolder).SOTA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SOFA(nsc) = 2 * Condition(ifolder).SOFA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SPPA(nsc) = 2 * Condition(ifolder).SPPA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SPTA(nsc) = 2 * Condition(ifolder).SPTA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SPFA(nsc) = 2 * Condition(ifolder).SPFA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).STTA(nsc) = 2 * Condition(ifolder).STTA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).STFA(nsc) = 2 * Condition(ifolder).STFA(nsc) / (nsc*nsc);   % normalization
        Condition(ifolder).SFFA(nsc) = 2 * Condition(ifolder).SFFA(nsc) / (nsc*nsc);   % normalization
    end

end
    
run = 1;
plotTwoCumulative(Condition(1).SOOA,Condition(5).SOOA,'HighAttention','RestingState','O-O',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(1).SPPA,Condition(5).SPPA,'HighAttention','RestingState','P-P',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(1).SFFA,Condition(5).SFFA,'HighAttention','RestingState','F-F',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(1).SOPA,Condition(5).SOPA,'HighAttention','RestingState','O-P',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(1).SPFA,Condition(5).SPFA,'HighAttention','RestingState','P-F',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(1).SOFA,Condition(5).SOFA,'HighAttention','RestingState','O-F',NA,run,settings,preprocessed_step);

run = 2;
plotTwoCumulative(Condition(2).SOOA,Condition(6).SOOA,'HighAttention','RestingState','O-O',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(2).SPPA,Condition(6).SPPA,'HighAttention','RestingState','P-P',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(2).SFFA,Condition(6).SFFA,'HighAttention','RestingState','F-F',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(2).SOPA,Condition(6).SOPA,'HighAttention','RestingState','O-P',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(2).SPFA,Condition(6).SPFA,'HighAttention','RestingState','P-F',NA,run,settings,preprocessed_step);

plotTwoCumulative(Condition(2).SOFA,Condition(6).SOFA,'HighAttention','RestingState','O-F',NA,run,settings,preprocessed_step);

  
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

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Four-lobes','-',label1,'-',label2,'-',titlelabel,'-','Run','-',int2str(run),'-',preprocessed_step,'.pdf'));
   

end
