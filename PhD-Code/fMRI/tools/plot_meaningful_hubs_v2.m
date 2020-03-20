
clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

Hubs = [264 266 511 358 573];

nClusters = 758;
nROI = 90;

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

marker{1} = 'b.';
marker{2} = 'b.';
marker{3} = 'r<';
marker{4} = 'r<';
marker{5} = 'k>';
marker{6} = 'k>';

line_color(1,:) = [0 0 255]./255;
line_color(2,:) = [0 0 255]./255;
line_color(3,:) = [255 0 0]./255;
line_color(4,:) = [255 0 0]./255;
line_color(5,:) = [0 0 0]./255;
line_color(6,:) = [0 0 0]./255;

line_style{1} = '-';
line_style{2} = '--';
line_style{3} = '-';
line_style{4} = '--';
line_style{5} = '-';
line_style{6} = '--';

marker_size(1) = 10;
marker_size(2) = 10;
marker_size(3) = 10;
marker_size(4) = 10;
marker_size(5) = 10;
marker_size(6) = 10;

% FC_strength = zeros(1,length(MeaningFul_Attention)-1);
% GC_FROM_strength = zeros(1,length(MeaningFul_Attention)-1);
% GC_TO_strength = zeros(1,length(MeaningFul_Attention)-1); 
% 
% for iMean=2:length(MeaningFul_Attention)
%     
%     if ~isempty(length(MeaningFul_Attention{iMean,7}))
%         FC_strength(iMean-1) = length(MeaningFul_Attention{iMean,7}); % FC - direct
%     else
%         FC_strength(iMean-1) = 0;
%     end
%     
%     if ~isempty(length(MeaningFul_Attention{iMean,8}))
%         FC_strength(iMean-1) = FC_strength(iMean-1) + length(MeaningFul_Attention{iMean,8}); % FC - indirect
%     else
%         FC_strength(iMean-1) = FC_strength(iMean-1) + 0;
%     end
% 
%     if ~isempty(length(MeaningFul_Attention{iMean,11}))
%         GC_FROM_strength(iMean-1) = length(MeaningFul_Attention{iMean,11}); % GC - direct - FROM
%     else
%         GC_FROM_strength(iMean-1) = 0;
%     end
%     
%     if ~isempty(length(MeaningFul_Attention{iMean,12}))
%         GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + length(MeaningFul_Attention{iMean,12}); % GC - indirect - FROM
%     else
%         GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + 0;
%     end
%     
%     if ~isempty(length(MeaningFul_Attention{iMean,15}))
%         GC_TO_strength(iMean-1) = length(MeaningFul_Attention{iMean,15}); % GC - direct - TO
%     else
%         GC_TO_strength(iMean-1) = 0;
%     end
%     
%     if ~isempty(length(MeaningFul_Attention{iMean,16}))
%         GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + length(MeaningFul_Attention{iMean,16}); % GC - indirect - TO
%     else
%         GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + 0;
%     end
%     
% end
% 
% [s i] = sort(FC_strength,'descend');
% idx_cluster(1) = i(1)+1;
% 
% [s i] = sort(GC_FROM_strength,'descend');
% idx_cluster(2) = i(1)+1;
% 
% [s i] = sort(GC_TO_strength,'descend');
% idx_cluster(3) = i(1)+1;
% 
% idx_cluster

%%% STRONGEST HUBS: 266 (FC, GC_TO) , 511 (GC_FROM)

for iHub=Hubs
% for iHub=264
    
    con_leng{1} = 'FC - direct';
    con_leng{2} = 'FC - indirect';
    con_leng{3} = 'GC-FROM - direct';
    con_leng{4} = 'GC-FROM - indirect';
    con_leng{5} = 'GC-TO - direct';
    con_leng{6} = 'GC-TO - indirect';
    
    all_clusters = [];
   
    for iMean=2:length(MeaningFul_Attention)
    
        if MeaningFul_Attention{iMean,1} == iHub
            
            hub_label = MeaningFul_Attention{iMean,2};
           
            conn_clusters(1).con = MeaningFul_Attention{iMean,7}; % FC - direct
            
            conn_clusters(2).con = MeaningFul_Attention{iMean,8}; % FC - indirect
            
            conn_clusters(3).con = MeaningFul_Attention{iMean,11}; % GC - direct - FROM
            
            conn_clusters(4).con = MeaningFul_Attention{iMean,12}; % GC - indirect - FROM
            
            conn_clusters(5).con = MeaningFul_Attention{iMean,15}; % GC - direct - TO
            
            conn_clusters(6).con = MeaningFul_Attention{iMean,16}; % GC - indirect - TO
            
            all_clusters = [conn_clusters(1).con,conn_clusters(2).con,conn_clusters(3).con,conn_clusters(4).con,conn_clusters(5).con',conn_clusters(6).con'];
            
        end
        
    end
    
    f = figure;

    r = 1;
    angles = linspace(0,2*pi,nClusters+1);
    x = r.*cos(angles);
    y = r.*sin(angles);
    
    step = angles(2) - angles(1);

    x_up = (r+0.05).*cos(angles);
    y_up = (r+0.05).*sin(angles);

    x_down = (r-0.05).*cos(angles);
    y_down = (r-0.05).*sin(angles);

    %plot(x,y,'k');
    %hold on
    plot(x,y,'k.','MarkerSize',3);
    hold on
    
    plot(x(iHub),y(iHub),'k.','MarkerSize',12);
    %plot(0,0,'k.','MarkerSize',30);

    for iAAL=1:nROI
    
        if iAAL == 1
            position = ceil(ROI_info{iAAL,5}/2);
        else
            position = ROI_info{iAAL,4} + round((ROI_info{iAAL,5} - ROI_info{iAAL,4})/2);
        end
        if mod(ROI_info{iAAL,5} - ROI_info{iAAL,4} + 1,2) ~= 0
            step = 0;
        else
            step = angles(2) - angles(1);
        end
        x_text = (r+0.20).*cos(angles+step/2);
        y_text = (r+0.20).*sin(angles+step/2);
        %h = text(x_text(position),y_text(position),ROI_info{iAAL,2},'FontSize',5);
        %set(h,'Rotation',(angles(position)+step/2)*180/pi');
        
        step = angles(2) - angles(1);
        x_div = r.*cos(angles - step/2);
        y_div = r.*sin(angles - step/2);
        x_div_up = (r+0.05).*cos(angles - step/2);
        y_div_up = (r+0.05).*sin(angles - step/2);
        plot([x_div(ROI_info{iAAL,4}) x_div_up(ROI_info{iAAL,4})],[y_div(ROI_info{iAAL,4}) y_div_up(ROI_info{iAAL,4})],'k');
        hold on
        
    end

    for iiCluster=1:nClusters

        idx_cluster = 0;
        seeds = [];
        for iROI=1:nROI

            for iCluster=1:length(ROI(iROI).clusters)

                idx_cluster = idx_cluster + 1;

                if idx_cluster == iiCluster

                    idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

                    Voxels = zeros(length(idx_voxels),3);

                    for iVoxel=1:length(idx_voxels)

                        [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

                    end

                    seeds = seeIfAVoxelIsInsideASeed(Voxels);

                end

            end

        end

         if ~isempty(seeds)

            if ~iscell(seeds)

                seed_to_plot = seeds;

            else

                seed_to_plot = seeds{1};

                for iseed=2:length(seeds)

                    seed_to_plot = strcat(seed_to_plot,' / ',seeds{iseed});

                end

            end

            x_seed = (r+0.05).*cos(angles(iiCluster));
            y_seed = (r+0.05).*sin(angles(iiCluster));
            %h = text(x_seed,y_seed,seed_to_plot,'FontSize',2.5);
            %set(h,'Rotation',angles(iiCluster)*180/pi);

        end

    end

    
    for icon=1:6
        
        index_to_remove = [];
        
        if ~isempty(conn_clusters(icon).con)
        
            for iicon=1:length(conn_clusters(icon).con)

                plot(x(conn_clusters(icon).con(iicon)),y(conn_clusters(icon).con(iicon)),marker{icon},'MarkerSize',marker_size(icon));
                hold on
                line([x(conn_clusters(icon).con(iicon)) x(iHub)],[y(conn_clusters(icon).con(iicon)) y(iHub)],'Color',line_color(icon,:),'LineStyle',line_style{icon},'LineWidth',2);
                
            end

        else
            
            index_to_remove = [index_to_remove;icon];
            
        end
        
    end

    con_leng(index_to_remove) = [];
        
    legend(con_leng);
    
    axis off
    
    print(f,'-depsc',strcat('Hubs-Attention-',hub_label,'-',int2str(iHub),'without-labels.eps'));
    
end
    