
clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

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

nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

for iDTI=1:8
  
    C = DTI(iDTI).C; 
  
    f = figure;

    imagesc(C);

    title(strcat('Subject-',int2str(iDTI)));

    clrmp = colormap('jet');
    clrmp(1,:) = [1 1 1];

    min_C = 0;
    max_C = 5000;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);

    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('DTI-Subject-',int2str(iDTI)));
    
end
    