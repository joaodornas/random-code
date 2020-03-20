function lowhigh_FC_ICA_occipital_parietal_lobe


settings_jan_0805;
%settings_elena_2905;

preprocessed_step = 'residual';
%preprocessed_step = 'filtered';

%doTheMath(settings,preprocessed_step);

%plotResults(settings,preprocessed_step);

%plotCumulative(settings,preprocessed_step);

settings_jan_0805;
all_subjects(1).settings = settings;
clear settings;
settings_elena_2905;
all_subjects(2).settings = settings;
clear settings;
plotAllSubjectsCumulative(all_subjects,preprocessed_step);

end

function doTheMath(settings,preprocessed_step)

preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed);

%ICA_preprocessed = 'ICA_lobe';
ICA_preprocessed = 'ICA_lobe_filtered';

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
    
    ICA_n(ifolder).nOccipital = size(occipital_ic,2);
    ICA_n(ifolder).nParietal = size(parietal_ic,2);
    
    allComponents = [occipital_ic,parietal_ic];
    
    [ICA_FC(ifolder).rho,ICA_FC(ifolder).pval] = corr(allComponents);
  
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-lobe-All-Runs-',preprocessed_step,'.mat'),'ICA_label','ICA_folder','ICA_n','ICA_FC','ICA_preprocessed');

end


function plotResults(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-lobe-All-Runs-',preprocessed_step,'.mat'));

for ifolder=1:length(ICA_folder)
    
   max_rho(ifolder) = max(max(ICA_FC(ifolder).rho));
   min_rho(ifolder) = min(min(ICA_FC(ifolder).rho));
   
%    max_percent_rho(ifolder) =  max(max(ICA_FC(ifolder).percent_rho));
%    min_percent_rho(ifolder) = min(min(ICA_FC(ifolder).percent_rho));
    
end

max_all = max(max_rho);
min_all = min(min_rho);

% max_percent_all = min(max_percent_rho);
% min_percent_all = min(min_percent_rho);

for ifolder=1:length(ICA_folder)
    
   f = figure;
   
   imagesc(ICA_FC(ifolder).rho,[min_all max_all]);
   
   xlabel('occipital - parietal');
   ylabel('parietal - occipital');
   
   title(ICA_label{ifolder});
   
   colorbar;
   
   hold on;
   
   p1(1) = ICA_n(ifolder).nOccipital;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nOccipital;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nOccipital;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   p2(2) = ICA_n(ifolder).nOccipital;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.pdf'));
   
end

end


function plotCumulative(settings,preprocessed_step)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-lobe-All-Runs-',preprocessed_step,'.mat'));

for ifolder=1:6
    
    NoA = ICA_n(ifolder).nOccipital;
    NpA = ICA_n(ifolder).nParietal;
   
    OOA = ICA_FC(ifolder).rho(1:NoA,1:NoA);

    PPA = ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);

    OPA = ICA_FC(ifolder).rho(1:NoA, NoA+1:NoA+NpA); 
    
    NA = min( [NoA NpA] );  

    SOOA = zeros(1,NA);
    SPPA = zeros(1,NA);
    SOPA = zeros(1,NA);

    for nsc = 2:NA    % number of components to be summed

        for ic = 2 : nsc   % loop over components in i-direction
            for jc = 1 : ic-1   % loop over componnets in j-direction

                SOOA(nsc) = SOOA(nsc) + abs( OOA(ic,jc) );
                SPPA(nsc) = SPPA(nsc) + abs( PPA(ic,jc) );
                SOPA(nsc) = SOPA(nsc) + abs( OPA(ic,jc) );

            end
        end

        SOOA(nsc) = 2 * SOOA(nsc) / (nsc*nsc);   % normalization
        SPPA(nsc) = 2 * SPPA(nsc) / (nsc*nsc);   % normalization
        SOPA(nsc) = 2 * SOPA(nsc) / (nsc*nsc);   % normalization
        
    end
    
    f = figure;
    
    fs = 14;
    
    hold on;
    plot(1:NA, SOOA, 'k', 'LineWidth', 2 );
    plot(1:NA, SPPA, 'b', 'LineWidth', 2 );
    plot(1:NA, SOPA, 'r', 'LineWidth', 2 );
    hold off;
    
    h = legend('O-O', 'P-P', 'O-P');
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(ICA_label{ifolder}, 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'-',preprocessed_step,'.pdf'));
   
end

end

function plotAllSubjectsCumulative(all_subjects,preprocessed_step)

experiment = all_subjects(1).settings.folders.experiment;

all_NA = [];
for isubject=1:length(all_subjects)
    
    settings = all_subjects(isubject).settings;

    analysis = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-lobe-All-Runs-',preprocessed_step,'.mat'));

    for ifolder=1:6
        
        NoA = analysis.ICA_n(ifolder).nOccipital;
        NpA = analysis.ICA_n(ifolder).nParietal;
        
        all_NA = [all_NA, NoA, NpA];
    
    end
    
    clear settings;
    clear analysis;
    
end

NA = min(all_NA);
    
for isubject=1:length(all_subjects)
    
    settings = all_subjects(isubject).settings;

    analysis = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-lobe-All-Runs-',preprocessed_step,'.mat'));

    for ifolder=1:6
  
        NoA = analysis.ICA_n(ifolder).nOccipital;
        NpA = analysis.ICA_n(ifolder).nParietal;
   
        OOA = analysis.ICA_FC(ifolder).rho(1:NoA,1:NoA);

        PPA = analysis.ICA_FC(ifolder).rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);

        OPA = analysis.ICA_FC(ifolder).rho(1:NoA, NoA+1:NoA+NpA); 
    
        SOOA = zeros(1,NA);
        SPPA = zeros(1,NA);
        SOPA = zeros(1,NA);

        for nsc = 2:NA    % number of components to be summed

            for ic = 2 : nsc   % loop over components in i-direction
                for jc = 1 : ic-1   % loop over componnets in j-direction

                    SOOA(nsc) = SOOA(nsc) + abs( OOA(ic,jc) );
                    SPPA(nsc) = SPPA(nsc) + abs( PPA(ic,jc) );
                    SOPA(nsc) = SOPA(nsc) + abs( OPA(ic,jc) );

                end
            end

            SOOA(nsc) = 2 * SOOA(nsc) / (nsc*nsc);   % normalization
            SPPA(nsc) = 2 * SPPA(nsc) / (nsc*nsc);   % normalization
            SOPA(nsc) = 2 * SOPA(nsc) / (nsc*nsc);   % normalization

        end
        
        subjects(isubject).folder(ifolder).SOOA = SOOA;
        subjects(isubject).folder(ifolder).SPPA = SPPA;
        subjects(isubject).folder(ifolder).SOPA = SOPA;
    
    end
    
    clear settings;
    clear analysis;
    
end

High_SOOA = zeros(1,NA);
High_SPPA = zeros(1,NA);
High_SOPA = zeros(1,NA);

Low_SOOA = zeros(1,NA);
Low_SPPA = zeros(1,NA);
Low_SOPA = zeros(1,NA);

Rest_SOOA = zeros(1,NA);
Rest_SPPA = zeros(1,NA);
Rest_SOPA = zeros(1,NA);

for isubject=1:length(all_subjects)
    
    High_SOOA = High_SOOA + subjects(isubject).folder(1).SOOA + subjects(isubject).folder(2).SOOA;
    Low_SOOA = Low_SOOA + subjects(isubject).folder(3).SOOA + subjects(isubject).folder(4).SOOA;
    Rest_SOOA = Rest_SOOA + subjects(isubject).folder(5).SOOA + subjects(isubject).folder(6).SOOA;

    High_SPPA = High_SPPA + subjects(isubject).folder(1).SPPA + subjects(isubject).folder(2).SPPA;
    Low_SPPA = Low_SPPA + subjects(isubject).folder(3).SPPA + subjects(isubject).folder(4).SPPA;
    Rest_SPPA = Rest_SPPA + subjects(isubject).folder(5).SPPA + subjects(isubject).folder(6).SPPA;

    High_SOPA = High_SOPA + subjects(isubject).folder(1).SOPA + subjects(isubject).folder(2).SOPA;
    Low_SOPA = Low_SOPA + subjects(isubject).folder(3).SOPA + subjects(isubject).folder(4).SOPA;
    Rest_SOPA = Rest_SOPA + subjects(isubject).folder(5).SOPA + subjects(isubject).folder(6).SOPA;

end

High_SOOA = High_SOOA ./ (2*length(all_subjects));
High_SPPA = High_SPPA ./ (2*length(all_subjects));
High_SOPA = High_SOPA ./ (2*length(all_subjects));

Low_SOOA = Low_SOOA ./ (2*length(all_subjects));
Low_SPPA = Low_SPPA ./ (2*length(all_subjects));
Low_SOPA = Low_SOPA ./ (2*length(all_subjects));

Rest_SOOA = Rest_SOOA ./ (2*length(all_subjects));
Rest_SPPA = Rest_SPPA ./ (2*length(all_subjects));
Rest_SOPA = Rest_SOPA ./ (2*length(all_subjects));

f = figure;

fs = 14;

hold on;
plot(1:NA, High_SOOA, 'k', 'LineWidth', 2 );
plot(1:NA, High_SPPA, 'b', 'LineWidth', 2 );
plot(1:NA, High_SOPA, 'r', 'LineWidth', 2 );
hold off;

h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
axis 'square';
axis([0 NA 0 0.5]);
title(strcat('High Attention','-',int2str(length(all_subjects)),'-','subjects'), 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );

print(f,'-djpeg',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','High-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.jpeg'));
print(f,'-depsc',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','High-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.eps'));
print(f,'-dpdf',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','High-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.pdf'));

f = figure;

fs = 14;

hold on;
plot(1:NA, Low_SOOA, 'k', 'LineWidth', 2 );
plot(1:NA, Low_SPPA, 'b', 'LineWidth', 2 );
plot(1:NA, Low_SOPA, 'r', 'LineWidth', 2 );
hold off;

h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
axis 'square';
axis([0 NA 0 0.5]);
title(strcat('Low Attention','-',int2str(length(all_subjects)),'-','subjects'), 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );

print(f,'-djpeg',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','Low-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.jpeg'));
print(f,'-depsc',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','Low-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.eps'));
print(f,'-dpdf',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','Low-Attention','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.pdf'));

f = figure;

fs = 14;

hold on;
plot(1:NA, Rest_SOOA, 'k', 'LineWidth', 2 );
plot(1:NA, Rest_SPPA, 'b', 'LineWidth', 2 );
plot(1:NA, Rest_SOPA, 'r', 'LineWidth', 2 );
hold off;

h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
axis 'square';
axis([0 NA 0 0.5]);
title(strcat('Resting State','-',int2str(length(all_subjects)),'-','subjects'), 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );

print(f,'-djpeg',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','RestingState','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.jpeg'));
print(f,'-depsc',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','RestingState','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.eps'));
print(f,'-dpdf',strcat(experiment,'-','Cumulative-ICA-Occipital-Parietal-lobe','-','RestingState','-',int2str(length(all_subjects)),'-','subjects','-',preprocessed_step,'.pdf'));

end
  
