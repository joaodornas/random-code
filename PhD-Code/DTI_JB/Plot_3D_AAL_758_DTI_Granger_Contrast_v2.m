
clear all

% parameters
nClusters = 758;
nROI = 90;
zcriterion = 2.3;
% zcriterion = 3.0;
nConn = 4;

% load AAL
load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:2:nROI-1
    
    roi_label{iROI+1} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
end

for iROI=2:2:nROI
    
    roi_label{iROI-1} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
end

% load 758 Clusters Info
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

area(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % frontal
area(2).idx = [27 28 29 30 31 32 33 34 35 36]; % rec_ins_cing
area(3).idx = [37 38 39 40 41 42]; % hc_amyg
area(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % occipital
area(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % parietal
area(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % temporal
area(7).idx = [71 72 73 74 75 76 77 78]; % subcortical

area(1).label = 'frontal';
area(2).label = 'rec_ins_cing';
area(3).label = 'hc_amyg';
area(4).label = 'occipital';
area(5).label = 'parietal';
area(6).label = 'temporal';
area(7).label = 'subcortical';

area(1).color = 'r';
area(2).color = 'b';
area(3).color = 'g';
area(4).color = 'y';
area(5).color = 'w';
area(6).color = 'm';
area(7).color = 'k';

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

        
for iConn=1:nConn

    C = Connectome(iConn).C;

    C_new = zeros(size(C));

    for iCluster=1:nClusters

        connections = squeeze(C(iCluster,:));

        abs_connections = abs(connections);

        [s,i] = sort(abs_connections,'descend');

        C_new(iCluster,i(1)) = connections(i(1));

    end           
             
    for iROI=1:nROI
   
        causal_from = C_new(ROI_info{iROI,4}:ROI_info{iROI,5},:);
        
        [s_from,i_from] = sort(abs(causal_from(:)),'descend');
     
        max_all = causal_from(i_from(1));
        [x, iCluster] = ind2sub(size(C_new(ROI_info{iROI,4}:ROI_info{iROI,5},:)),i_from(1));
      
        for iiiROI=1:nROI
            
            if iCluster >= ROI_info{iiiROI,4} && iCluster <= ROI_info{iiiROI,5}; iiROI = iiiROI; end
            
        end

        Connectome(iConn).C_AAL(iROI,iiROI) = max_all;

    end
    
end

for iConn=[1 3]
    
    C = Connectome(iConn).C_AAL;

    C = C./max(abs(C(:)));

    load aal_cog.txt % Only 90 areas (not 116)
    MNI_coord=aal_cog;
    clear aal_cog

    N=length(C);

    figure('name',Connectome(iConn).label)
    set(gcf, 'color', 'white');

    hold on

    for n=1:N

        for p = 1:N

            if C(n,p) > 0; color_label = 'r'; elseif C(n,p) < 0; color_label = 'b'; end

            direction = MNI_coord(p,:) - MNI_coord(n,:);

            if C(n,p) ~= 0

                quiver3(MNI_coord(n,1),MNI_coord(n,2),MNI_coord(n,3),direction(1),direction(2),direction(3),color_label,'LineWidth',2);

            end

        end

    end

    % PLOT THE NODES AS 3D SPHERES
    [x,y,z] = sphere;
    sphere_size = 2;
    for n=1:N

        for iArea=1:7

            if ismember(n,area(iArea).idx); color_label = area(iArea).color; end

        end
        
        x_text = MNI_coord(n,1);
        y_text = MNI_coord(n,2);
        
        if n == 26 || n == 8 || n == 43 || n == 83; y_text = y_text + 1; end
        if n == 13 || n == 64 ; y_text = y_text - 1; end
        if n == 75 || n == 76; y_text = y_text - 5; end
        if n == 12; y_text = y_text + 2; end
        if n == 84 || n == 87; y_text = y_text + 5; end
        if n == 27 || n == 29 || n == 30 || n == 73 || n == 74 || n == 85 || n == 90; y_text = y_text + 3; end
        if n == 63 || n == 77 || n == 78; y_text = y_text - 2; end
        if n == 28 || n == 44 || n == 50; y_text = y_text -5; end
        if n == 31 || n == 70 || n == 71 || n == 88; y_text = y_text - 3; end
        
        if n == 26 || n == 35 || n == 36; x_text = x_text + 5; end
        
        surf(sphere_size*x+MNI_coord(n,1),sphere_size*y+MNI_coord(n,2),sphere_size*z+MNI_coord(n,3),'FaceColor',color_label,'EdgeColor','none');
        text(x_text,y_text,sphere_size*max(z(:))+MNI_coord(n,3),roi_label{n},'FontSize',8);
        
    end
    axis off;
    axis equal

    camlight;
    rotate3d;
    material dull; 
    lighting phong;

    %----------------- To plot the cortex/scalp underneath-------------------------

    % put your own spm path 
    addpath 'Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8'
    addpath 'Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical'
    D=spm_eeg_load;
    % Only forthis template
    fn = fieldnames(D.inv{1}.mesh);
    tess_ind = find(strncmp(fn,'tess_',5) & ~strcmp(fn,'tess_mni'));
    for i = 1:length(tess_ind)
        location=D.inv{1}.mesh.(fn{tess_ind(i)});
        D.inv{1}.mesh.(fn{tess_ind(i)})=location(14:end);
    end
    mesh = spm_eeg_inv_transform_mesh(D.inv{1}.mesh.Affine, D.inv{1}.mesh);
    Mscalp  = mesh.tess_scalp;
    headface    = Mscalp.face;
    headvert    = Mscalp.vert;
    cortexface  = mesh.tess_ctx.face;
    cortexvert  = mesh.tess_ctx.vert;

    %plot3(headface(:,1),headface(:,2),headface(:,3),'.','MarkerSize',1,'Color',[0.95 .7 .55])
    %patch('vertices',headvert,'faces',headface,'EdgeColor','none','FaceColor',[1 .7 .55],'FaceAlpha',0.1); %plots a mesh surface for the headshape
    patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.5 0.5 .5],'FaceAlpha',0.15); %plots a mesh surface for the cortex

    camlight;
    rotate3d;
    material dull; 
    lighting phong;
    
end