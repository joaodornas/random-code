
% parameters
nClusters = 758;
nROI = 90;
zcriterion = 2.3;
nConn = 4;

% load 758 Clusters Info
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

area(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % frontal
area(2).idx = [27 28 29 30 31 32 33 34 35 36]; % rec_ins_cing
area(3).idx = [37 38 39 40 41 42]; % hc_amyg
area(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % occipital
area(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % parietal
area(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % temporal
area(7).idx = [71 72 73 74 75 76 77 78]; % subcortical

for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

% load Granger Contrast
load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));

% compute Attention & Stimulus Only Results
Attention_rho = Attention_Contrast.Z;
Attention_rho(find(Attention_rho>(-1)*zcriterion & Attention_rho<zcriterion)) = 0; 

Stimulus_rho = Stimulus_Contrast.Z;
Stimulus_rho(find(Stimulus_rho>(-1)*zcriterion & Stimulus_rho<zcriterion)) = 0; 

Attention_rho(isnan(Attention_rho)) = 0;
Stimulus_rho(isnan(Stimulus_rho)) = 0;

AttentionOnly_rho = Attention_rho; 
AttentionOnly_rho(find(Stimulus_rho)) = 0;

StimulusOnly_rho = Stimulus_rho; 
StimulusOnly_rho(find(Attention_rho)) = 0; 

% load DTI
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

% compute Common DTI
ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

Connectome(1).C = AttentionOnly_rho;
Connectome(1).C(~(common_DTI)) = 0;
Connectome(1).label = 'Attention-Direct';
Connectome(2).C = AttentionOnly_rho;
Connectome(2).C(find(common_DTI)) = 0;
Connectome(2).label = 'Attention-Indirect';
Connectome(3).C = StimulusOnly_rho;
Connectome(3).C(~(common_DTI)) = 0;
Connectome(3).label = 'Stimulus-Direct';
Connectome(4).C = StimulusOnly_rho;
Connectome(4).C(find(common_DTI)) = 0;
Connectome(4).label = 'Stimulus-Indirect';

% load AAL
load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

rgb(1) = 192;
rgb(2) = 128;
rgb(3) = 32;
rgb(4) = 224;
rgb(5) = 160;
rgb = rgb./255;

for iROI=1:nROI
    
    roi_label{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    if ismember(iROI,area(1).idx)
        
        aal_color = 1;
        
    elseif ismember(iROI,area(2).idx)
        
        aal_color = 2;
        
    elseif ismember(iROI,area(3).idx)
        
        aal_color = 3;
        
    elseif ismember(iROI,area(4).idx)
        
        aal_color = 4;
        
    elseif ismember(iROI,area(5).idx)
        
        aal_color = 5;
        
     elseif ismember(iROI,area(6).idx)
        
        aal_color = 3;
        
     elseif ismember(iROI,area(7).idx)
        
        aal_color = 2;
        
    end
    
    roi_color(iROI,:) = [rgb(aal_color), rgb(aal_color), rgb(aal_color)];
    
end

theta = linspace(0,2*pi,nClusters+1);

r = 2;
x = r .* cos(theta);
y = r .* sin(theta);
r_max = 0;

x_AAL_labels = x .* ( r_max + 1);
y_AAL_labels = y .* ( r_max + 1);

for iROI=1:nROI
    
    idx_cluster(iROI) = round( ROI_info{iROI,4} + (ROI_info{iROI,5} - ROI_info{iROI,4})/2 );
    
    h = text(x_AAL_labels(idx_cluster(iROI))*1.01,y_AAL_labels(idx_cluster(iROI))*1.01,roi_label{iROI},'FontSize',5);
    set(h,'Rotation',(theta(idx_cluster(iROI))*180/pi));
    
end

hold on

plot(x,y,'k');
r = 1.8;
x = r .* cos(theta);
y = r .* sin(theta);
plot(x,y,'k');

xlim([-2.5 2.5]);
ylim([-2.5 2.5]);

