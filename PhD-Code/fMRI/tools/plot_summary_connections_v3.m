
clear all

% for parcellation_label = {'MD758', 'SP758'}
for parcellation_label = {'MD758'}

    for direction_label = {'Increase', 'Decrease'}

        load(strcat('FC-Voxels-',parcellation_label{1},'-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only-',direction_label{1},'-corrected.mat'));
    
        for contrast_label = {'Attention', 'Stimulus'}
            
            eval(strcat('MeaningFul_Clusters = MeaningFul_',contrast_label{1},';'));
            
                load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');
                
                clearvars -except ROI MeaningFul_Clusters contrast_label direction_label parcellation_label MeaningFul_Stimulus MeaningFul_Attention

                nClusters = 758;
                
                disp(strcat(parcellation_label{1},'-',direction_label{1},'-',contrast_label{1}));

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

                FC_strength = zeros(1,length(MeaningFul_Clusters)-1);
                GC_FROM_strength = zeros(1,length(MeaningFul_Clusters)-1);
                GC_TO_strength = zeros(1,length(MeaningFul_Clusters)-1);

                do_indirect = 1;

                for iMean=2:length(MeaningFul_Clusters)

                    idx_cluster(iMean-1) = MeaningFul_Clusters{iMean,1};

                    cluster_label{iMean-1} = MeaningFul_Clusters{iMean,2};

                    if ~isempty(length(MeaningFul_Clusters{iMean,7}))
                        FC_strength(iMean-1) = length(MeaningFul_Clusters{iMean,7}); % FC - direct
                    else
                        FC_strength(iMean-1) = 0;
                    end
                    if ~isempty(length(MeaningFul_Clusters{iMean,8})) && do_indirect
                        FC_strength(iMean-1) = FC_strength(iMean-1) + length(MeaningFul_Clusters{iMean,8}); % FC - indirect
                    else
                        FC_strength(iMean-1) = FC_strength(iMean-1) + 0;
                    end

                    if ~isempty(length(MeaningFul_Clusters{iMean,11}))
                        GC_FROM_strength(iMean-1) = length(MeaningFul_Clusters{iMean,11}); % GC - direct - FROM
                    else
                        GC_FROM_strength(iMean-1) = 0;
                    end
                    if ~isempty(length(MeaningFul_Clusters{iMean,12})) && do_indirect
                        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + length(MeaningFul_Clusters{iMean,12}); % GC - indirect - FROM
                    else
                        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + 0;
                    end

                    if ~isempty(length(MeaningFul_Clusters{iMean,15}))
                        GC_TO_strength(iMean-1) = length(MeaningFul_Clusters{iMean,15}); % GC - direct - TO
                    else
                        GC_TO_strength(iMean-1) = 0;
                    end
                    if~isempty(length(MeaningFul_Clusters{iMean,16})) && do_indirect
                        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + length(MeaningFul_Clusters{iMean,16}); % GC - indirect - TO
                    else
                        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + 0;
                    end

                end

                FC_strength = ( FC_strength ./ max(FC_strength) ) * 20;

                GC_FROM_percent = GC_FROM_strength ./ (GC_FROM_strength + GC_TO_strength);
                GC_TO_percent = GC_TO_strength ./ (GC_FROM_strength + GC_TO_strength);
               
                [u_cluster_label, c_clusters] = count_unique(cluster_label);
                
                disp(strcat('nClusters:',int2str(size(MeaningFul_Clusters,1))));
                disp(strcat('GC-FROM:',num2str(length(find(GC_FROM_percent>0.5)))));
                disp(strcat('GC-TO:',num2str(length(find(GC_TO_percent>0.5)))));
                disp(strcat('GC-EQUAL:',num2str(length(find(GC_FROM_percent==0.5)))));
                disp(strcat('nROIs:',int2str(length(u_cluster_label))));

                % for iUnique=1:length(u_cluster_label)
                % 
                %     for iROI=1:nROI
                %            
                %        if strcmp(ROI(iROI).label,u_cluster_label(iUnique))
                %            
                %            n_clusters(iUnique) = ROI(iROI).nClusters;
                %        
                %        end
                %        
                %     end
                %     
                % end

                f = figure;

                r = 1;
                % angles = linspace(0,2*pi,sum(n_clusters)+1);
                angles = linspace(0,2*pi,sum(c_clusters)+1);
                angles(end) = [];
                x = r.*cos(angles);
                y = r.*sin(angles);

                plot(x,y,'w');
                hold on

                plot(x,y,'w*','MarkerSize',2);
                hold on

                % seed_size = 0.5;
                % for iMean=2:length(MeaningFul_Clusters)
                %     
                %     seeds = MeaningFul_Clusters{iMean,17};
                %     
                %     if ~isempty(seeds)
                %         
                %         if ~iscell(seeds)
                %             
                %             seed_to_plot = seeds;
                %             
                %         else
                %             
                %             seed_to_plot = seeds{1};
                %             
                %             for iseed=2:length(seeds)
                %                 
                %                 seed_to_plot = strcat(seed_to_plot,' / ',seeds{iseed});
                %                 
                %             end
                %             
                %         end
                %         
                %         x_seed = (r-0.40).*cos(angles(iMean-1));
                %         y_seed = (r-0.40).*sin(angles(iMean-1));
                %         %h = text(x_seed,y_seed,seed_to_plot,'FontSize',seed_size);
                %         %set(h,'Rotation',angles(iMean-1)*180/pi);
                %         
                %     end
                %     
                % end

                AAL_label_size = 5;

                for iUnique=1:length(u_cluster_label)

                    if iUnique > 1

                        % last_clusters = sum(n_clusters(1:iUnique-1));
                        last_clusters = sum(c_clusters(1:iUnique-1));

                    else

                        last_clusters = 0;

                    end

                    % position = round(n_clusters(iUnique)/2) + last_clusters;
                    position = round(c_clusters(iUnique)/2) + last_clusters;

                    x_label = (r+0.40).*cos(angles(position));
                    y_label = (r+0.40).*sin(angles(position));
% % %                     
% % %                     if ( cos(angles(position)) >= 0 && sin(angles(position)) >=0 ) || ( cos(angles(position)) >= 0 && sin(angles(position)) <=0 ) 
% % %                     
% % %                         h = text(x_label,y_label,u_cluster_label(iUnique),'FontSize',AAL_label_size);
% % %                         set(h,'Rotation',angles(position)*180/pi);
% % % 
% % %                     else
% % %                         
% % %                         h = text(x_label,y_label,u_cluster_label(iUnique),'FontSize',AAL_label_size);
% % %                         set(h,'Rotation',angles(position)*180/pi);
% % %                         
% % %                     end
                    
                    % position = n_clusters(iUnique) + last_clusters;
                    position = c_clusters(iUnique) + last_clusters;

                    step = angles(2) - angles(1);

                    x_division_up = (r+0.05).*cos(angles(position) + step/2);
                    x_division = r.*cos(angles(position) + step/2);
                    y_division_up = (r+0.05).*sin(angles(position) + step/2);
                    y_division = r.*sin(angles(position) + step/2);

                    plot([x_division x_division_up],[y_division y_division_up],'w','LineWidth',2);
                    hold on

                end

                % for iMean=1:length(idx_cluster)
                %     
                %     for iROI=1:nROI
                %         
                %         if idx_cluster(iMean) >= ROI_info{iROI,4} && idx_cluster(iMean) <= ROI_info{iROI,5}
                %         
                %             order_clusters(iMean) = idx_cluster(iMean) - ROI_info{iROI,4} + 1;
                %             
                %         end
                %         
                %     end
                %     
                % end

                reduction = 0.15;
                FC_width = 0.1;
                GC_percent = 1;
                GC_ball = 12;

                for iUnique=1:length(u_cluster_label)

                   for iCluster=1:c_clusters(iUnique)

                       if iUnique > 1

                           start = sum(c_clusters(1:iUnique-1));

                       else

                           start = 0;

                       end

                       % order = order_clusters(start+iCluster);

                       position_hubs = start+iCluster;
                       position = start+iCluster;

                %        if iUnique > 1
                %            
                %            % position = sum(n_clusters(1:iUnique-1)) + order;
                %            position = sum(c_clusters(1:iUnique-1));
                %        
                %        else
                %            
                %            position = iCluster;
                %            
                %        end

% % %                         x_up = (r+GC_FROM_percent(position_hubs)*reduction).*cos(angles(position));
% % %                         y_up = (r+GC_FROM_percent(position_hubs)*reduction).*sin(angles(position));
% % % 
% % %                         x_down = (r-GC_TO_percent(position_hubs)*reduction).*cos(angles(position));
% % %                         y_down = (r-GC_TO_percent(position_hubs)*reduction).*sin(angles(position));
                        
                        x_up = (r+FC_strength(position_hubs)*FC_width*reduction).*cos(angles(position));
                        y_up = (r+FC_strength(position_hubs)*FC_width*reduction).*sin(angles(position));

% % %                         x_down = (r-1*reduction).*cos(angles(position));
% % %                         y_down = (r-1*reduction).*sin(angles(position));

                        x_mid = (r).*cos(angles(position));
                        y_mid = (r).*sin(angles(position));

% % %                         if FC_strength(position_hubs) ~=0
% % %                             plot([x_up x_mid],[y_up y_mid],'k','LineWidth',FC_strength(position_hubs)*FC_width);
% % %                             hold on
% % %                         end

% % %                         if GC_FROM_percent(position_hubs) > GC_TO_percent(position_hubs)
% % %                             color_GC = 'r.';
% % %                             percent = GC_FROM_percent(position_hubs);
% % %                             plot(x_up,y_up,color_GC,'MarkerSize',GC_ball);
% % %                         elseif GC_FROM_percent(position_hubs) < GC_TO_percent(position_hubs)
% % %                             color_GC = 'g.';
% % %                             percent = GC_TO_percent(position_hubs);
% % %                             plot(x_down,y_down,color_GC,'MarkerSize',GC_ball);
% % %                         elseif GC_FROM_percent(position_hubs) == GC_TO_percent(position_hubs)
% % %                             color_GC = 'y.';
% % %                             percent = GC_TO_percent(position_hubs);
% % %                             plot(x_mid,y_mid,color_GC,'MarkerSize',GC_ball);
% % %                         end



                %         x_up = (r+GC_FROM_percent(position_hubs)*reduction+0.01).*cos(angles(position));
                %         y_up = (r+GC_FROM_percent(position_hubs)*reduction+0.01).*sin(angles(position));
                %         h = text(x_up,y_up,int2str(round(percent*100)),'FontSize',GC_percent);
                %         set(h,'Rotation',angles(position_hubs)*180/pi);

                   end

                end

                all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};
                color(1,:) = [0,0,128]./255;
                color(2,:) = [0.5,0,0.9];
                color(3,:) = [0,191,255]./255;
                color(4,:) = [0,100,0]./255;
                color(5,:) = [255,255,0]./255;
                color(6,:) = [0.91,0.41,0.17];
                color(8,:) = [255,0,0]./255;
                color(7,:) = [125,125,125]./255;
                
                for iROI=1:nROI
                    
                    ROI_label = ROI(iROI).label;
                    
                    this_ROI_is_here_cells = strfind(u_cluster_label,ROI_label);
                    
                    idx_ROI_inside_this_analysis = [];
                    value = [];
                    for i=1:length(this_ROI_is_here_cells)
                        value = cell2mat(this_ROI_is_here_cells(i));
                        if value
                            idx_ROI_inside_this_analysis = i;
                        end
                    end
                    
                    if idx_ROI_inside_this_analysis
                    
                        all_ROI_clusters = MeaningFul_Clusters(:,2);
                
                        all_clusters = strfind(all_ROI_clusters,ROI_label);
                
                        idx_clusters_this_ROI = [];
                        idx = 0;
                        for i=1:length(all_clusters)
                            if cell2mat(all_clusters(i))
                                idx = idx + 1;
                                idx_clusters_this_ROI(idx) = i;
                            end
                        end
                
                        seed_to_plot = [];
                        for iMean=idx_clusters_this_ROI
                
                            seeds = MeaningFul_Clusters{iMean,17};
                
                            if ~isempty(seeds)
                
                                if ~iscell(seeds)
                
                                    seed_to_plot = seeds;
                
                                else
                
                                    seed_to_plot = seeds{1};
                
                                    for iseed=2:length(seeds)
                
                                        seed_to_plot = strcat(seed_to_plot,' / ',seeds{iseed});
                
                                    end
                
                                end
                
                            end
                
                        end
                
                        if ~isempty(seed_to_plot)
                
                            idx_NET_to_plot = zeros(1,length(all_networks_labels));
                            for iNET=1:length(all_networks_labels)
                
                                NET_label = all_networks_labels{iNET};
                
                                comp = strfind(seed_to_plot,NET_label);
                                comp = length(find(comp));
                                
                                if comp
                                    idx_NET_to_plot(iNET) = 1;
                                end
                
                            end
                        
                            idx_all_NETs = find(idx_NET_to_plot);
                            
%                             raios = linspace(0.7,0.5,length(idx_all_NETs));
%                             raios = sort(raios);

                            raios = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9] + 0.2;
                            
                            for i=1:length(idx_all_NETs)
                
                                if idx_ROI_inside_this_analysis > 1
                                    
                                    start_cluster = sum(c_clusters(1:idx_ROI_inside_this_analysis-1));
                                    end_cluster = sum(c_clusters(1:idx_ROI_inside_this_analysis));
                                    
                                else
                                   
                                    start_cluster = 1;
                                    end_cluster = c_clusters(1);
                                    
                                end
                
                                x_patch_start_l = (r-raios(i)).*cos(angles(start_cluster));
                                x_patch_end_l = (r-raios(i)).*cos(angles(end_cluster));
                                y_patch_start_u = (r-raios(i)).*sin(angles(start_cluster));
                                y_patch_end_d = (r-raios(i)).*sin(angles(end_cluster));
                
                                x_patch_start_r = (r-raios(i)-0.1).*cos(angles(start_cluster));
                                x_patch_end_r = (r-raios(i)-0.1).*cos(angles(end_cluster));
                                y_patch_start_uu = (r-raios(i)-0.1).*sin(angles(start_cluster));
                                y_patch_end_dd = (r-raios(i)-0.1).*sin(angles(end_cluster));
                
%%%                                fill([x_patch_start_l x_patch_end_l x_patch_end_r x_patch_start_r],[y_patch_start_u y_patch_end_d y_patch_end_dd y_patch_start_uu],color(idx_all_NETs(i),:),'EdgeColor','none');
                                hold on
                
                            end
                
                        end
                    
                    end
                    
                end

                axis off

                print(f,'-depsc',strcat('Kiwis-',contrast_label{1},'-',direction_label{1},'-',parcellation_label{1},'-white-contour.eps'));
                
        end
        
    end
    
end




    