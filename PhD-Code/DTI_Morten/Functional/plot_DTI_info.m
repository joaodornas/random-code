
clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nSubjects = 8;

prefix = 'Z:\_DATA\LOW-HIGH-ATTENTION';
% sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX-SP758\Und_PRE_xx.mat';

% DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
% DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
% DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

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

ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:nSubjects
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ nSubjects;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

nConnections = 1000;

for iCon=nConnections:-1:1
    
    threshold = iCon;
    
    common_DTI_thresholded = common_DTI;
    common_DTI_thresholded(common_DTI_thresholded<threshold) = 0; 
    
    nC_Common(iCon) = length(find(common_DTI_thresholded));
    
    for iDTI=1:nSubjects
        
        DTI(iDTI).C_thresholded = DTI(iDTI).C;
        DTI(iDTI).C_thresholded(DTI(iDTI).C_thresholded<threshold) = 0;
        
    end
    
    DTI(1).unique = DTI(1).C_thresholded & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(2).unique = DTI(2).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(3).unique = DTI(3).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(4).unique = DTI(4).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(5).unique = DTI(5).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(6).unique = DTI(6).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(7).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(7).unique = DTI(7).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(8).C_thresholded);
    DTI(8).unique = DTI(8).C_thresholded & ~(DTI(1).C_thresholded) & ~(DTI(2).C_thresholded) & ~(DTI(3).C_thresholded) & ~(DTI(4).C_thresholded) & ~(DTI(5).C_thresholded) & ~(DTI(6).C_thresholded) & ~(DTI(7).C_thresholded);
               
    nC_Unique(iCon) = 0;
    
    for iDTI=1:8
        
       nC_Unique(iCon) = nC_Unique(iCon) + length(find(DTI(iDTI).unique)); 
        
    end
    
    any_DTI = DTI(1).C_thresholded | (DTI(2).C_thresholded) | (DTI(3).C_thresholded) | (DTI(4).C_thresholded) | (DTI(5).C_thresholded) | (DTI(6).C_thresholded) | (DTI(7).C_thresholded) | (DTI(8).C_thresholded);
    
    for iDTI=1:8
       
        any_DTI(find(DTI(iDTI).unique)) = 0;
        
    end
    
    any_DTI(find(common_DTI)) = 0;
    
    for iDTI=1:8
       
        DTI(iDTI).logical = 1 - ~(DTI(iDTI).C_thresholded);
        
    end
    
    weight_DTI = zeros(size(DTI(1).C));
    for iDTI=1:8
       
        weight_DTI = weight_DTI + DTI(iDTI).logical;
        
    end
    
    weight_DTI(find(~any_DTI)) = 0;
    
%     idx_any = find(any_DTI);
%     
%     for iAny=1:length(idx_any)
%         
%        weight(iAny) = 0;
%        
%        for iDTI=1:8
%        
%            if DTI(iDTI).C_thresholded(idx_any(iAny)) ~= 0; weight(iAny) = weight(iAny) + 1; end
%            
%        end
%        
%        %weight(iAny) = weight(iAny) ./ 8;
%        
%     end
    
    nC_Any(iCon) = length(find(any_DTI)) * mean(weight_DTI(:));
    
end

nC_Common = nC_Common * 8;

area_total = trapz(nC_Common);

for iC=1:length(nC_Common)
    
    area_parcial = trapz(nC_Common(1:iC));
    
    if round( (area_parcial / area_total) * 100 ) == 25; area_threshold(1) = iC; end
    if round( (area_parcial / area_total) * 100 ) == 50; area_threshold(2) = iC; end
    if round( (area_parcial / area_total) * 100 ) == 75; area_threshold(3) = iC; end
    
end


f = figure;
interval = 1000;
plot(nC_Common(1:interval),'b','LineWidth',2);
hold on

plot(nC_Unique(1:interval),'r','LineWidth',2);
plot(nC_Any(1:interval),'k','LineWidth',2);
legend({'Common', 'Unique', 'Any'});
xlabel('threshold');
ylabel('# of connections (absolute)');

h = area(1:area_threshold(1),nC_Common(1:area_threshold(1)));
set(h,'FaceColor',[0 0.5 1]);
g = area((area_threshold(1)+1):area_threshold(2),nC_Common((area_threshold(1)+1):area_threshold(2)));
set(g,'FaceColor',[1 0.5 0]);
w = area((area_threshold(2)):area_threshold(3),nC_Common((area_threshold(2)):area_threshold(3)));
set(w,'FaceColor',[1 0 0]);
k = area((area_threshold(3)+1):1000,nC_Common((area_threshold(3)+1):1000));
set(k,'FaceColor',[0 1 0]);

plot(nC_Unique(1:interval),'r','LineWidth',2);
plot(nC_Any(1:interval),'k','LineWidth',2);
legend({'Common', 'Unique', 'Any'});
xlabel('threshold');
ylabel('# of connections (absolute)');

plot(nC_Common(1:interval),'b','LineWidth',3);

for i=1:4
    
    if i == 1; start_thre = 20; end_thre = area_threshold(1); end
    if i == 2; start_thre = area_threshold(1); end_thre = area_threshold(2); end
    if i == 3; start_thre = area_threshold(2); end_thre = area_threshold(3); end
    if i == 4; start_thre = area_threshold(3); end_thre = 1000; end
    
    common_DTI_thresholded = common_DTI;
    common_DTI_thresholded(~(common_DTI_thresholded>start_thre & common_DTI_thresholded<=end_thre)) = 0; 

    common_DTI_thresholded(find(common_DTI_thresholded)) = 1;
    
    figure;
    imagesc(common_DTI_thresholded);
    
    if i == 1; cmap = colormap('parula'); cmap = flipud(cmap); max_v = 1; end
    if i == 2; cmap = colormap('jet'); max_v = 1.5; end
    if i == 3; cmap = colormap('jet'); max_v = 1.2; end
    if i == 4; cmap = colormap('winter'); max_v = 1; end
    
    cmap(1,:) = [1 1 1];
    caxis([0 max_v]);
    colormap(cmap);
    
end

overlaid = zeros(size(common_DTI_thresholded));
for i=1:4
    
    if i == 1; start_thre = 20; end_thre = area_threshold(1); color = 0.2; end
    if i == 2; start_thre = area_threshold(1); end_thre = area_threshold(2); color = 0.7; end
    if i == 3; start_thre = area_threshold(2); end_thre = area_threshold(3); color = 0.9; end
    if i == 4; start_thre = area_threshold(3); end_thre = 1000; color = 0.5; end
    
    common_DTI_thresholded = round(common_DTI);
    common_DTI_thresholded(~(common_DTI_thresholded>start_thre & common_DTI_thresholded<=end_thre)) = 0; 

    common_DTI_thresholded(find(common_DTI_thresholded)) = color;
    
    overlaid = overlaid + common_DTI_thresholded;
    
    imagesc(overlaid);
    caxis([0 1]);
    cmap = colormap('jet');
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    
end

%         jump = 0;
%         last_clusters_so_far = 0;
%         for iROI=1:nROI
% 
%            if ~isempty(find(ismember(x_idx,iROI)))
% 
%                jump = jump + 1;
% 
%                if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end
% 
%                if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end
% 
%                clusters_so_far = 0;
%                for iiROI=1:(limit_roi_clusters+roi_distance)
% 
%                    clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;
% 
%                end
% 
%                x_idx_jump(jump) = clusters_so_far;
% 
%            end
% 
%         end
% 
%         clusters_so_far = 0;
%         nClusters = 758;
%         for iROI=1:nROI
% 
%            clusters_so_far = ROI(iROI).nClusters + clusters_so_far;
% 
%            if ~isempty(find(ismember(x_idx,iROI)))
% 
%                hold on;
% 
%                plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');
% 
%                plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');
% 
%            end
% 
%         end
% 
%         ax = gca;
%         %ax.XTick = x_idx_jump;
%         %ax.XTickLabel = x_label;
%         set(ax,'XTick',x_idx_jump);
%         set(ax,'XTickLabel',x_label);
% 
%         %ax.YTick = x_idx_jump;
%         %ax.YTickLabel = x_label;
%         set(ax,'YTick',x_idx_jump);
%         set(ax,'YTickLabel',x_label);
% 
%         xticklabel_rotate;


