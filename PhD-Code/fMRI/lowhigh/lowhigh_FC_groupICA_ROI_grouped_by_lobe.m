function lowhigh_FC_groupICA_ROI_grouped_by_lobe


settings_jan_0805;
%settings_elena_2905;

%FC_groupICA_ROI_grouped_by_lobe(settings);
% FC_groupICA_ROI_grouped_by_lobe_Half(settings);

%FC_groupICA_ROI_grouped_by_lobe_clusters(settings);

%plot_FC_ICs(settings);

%plot_FC_SumOfSquare(settings);

% plot_FC_SumOfSquare_Together(settings);
% 
% clear settings
% settings_jan_0805;
% all_settings(1).settings = settings;
% clear settings
% settings_elena_2905;
% all_settings(2).settings = settings;
% clear settings

%plot_FC_lobes(all_settings);

%plot_FC_crosslobes(all_settings)

plot_FC_SumOfSquare_Together_Half(settings);



end

function FC_groupICA_ROI_grouped_by_lobe(settings)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

[FC_ICA_FO, FC_n_FO] = getFC_ICA(settings,idx_frontal,idx_occipital,'Frontal','Occipital');
[FC_ICA_PO, FC_n_PO] = getFC_ICA(settings,idx_parietal,idx_occipital,'Parietal','Occipital');
[FC_ICA_FP, FC_n_FP] = getFC_ICA(settings,idx_frontal,idx_parietal,'Frontal','Parietal');

[FC_SQ_FO] = getFC_SQ(FC_ICA_FO,FC_n_FO);
[FC_SQ_PO] = getFC_SQ(FC_ICA_PO,FC_n_PO);
[FC_SQ_FP] = getFC_SQ(FC_ICA_FP,FC_n_FP);

%save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'),'FC_ICA_FO','FC_n_FO');
save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'),'FC_ICA_FO','FC_n_FO','FC_ICA_PO','FC_n_PO','FC_ICA_FP','FC_n_FP','FC_SQ_FO','FC_SQ_PO','FC_SQ_FP');

end

function FC_groupICA_ROI_grouped_by_lobe_Half(settings)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

[FC_ICA_FO, FC_n_FO] = getFC_ICA_Half(settings,idx_frontal,idx_occipital,'Frontal','Occipital');
[FC_SQ_FO] = getFC_SQ_Half(FC_ICA_FO,FC_n_FO);
save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FO','.mat'),'FC_ICA_FO','FC_n_FO','FC_SQ_FO','-v7.3');

[FC_ICA_PO, FC_n_PO] = getFC_ICA_Half(settings,idx_parietal,idx_occipital,'Parietal','Occipital');
[FC_SQ_PO] = getFC_SQ_Half(FC_ICA_PO,FC_n_PO);
save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-PO','.mat'),'FC_ICA_PO','FC_n_PO','FC_SQ_PO','-v7.3');

[FC_ICA_FP, FC_n_FP] = getFC_ICA_Half(settings,idx_frontal,idx_parietal,'Frontal','Parietal');
[FC_SQ_FP] = getFC_SQ_Half(FC_ICA_FP,FC_n_FP);
save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FP','.mat'),'FC_ICA_FP','FC_n_FP','FC_SQ_FP','-v7.3');

%save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half','.mat'),'FC_ICA_FO','FC_n_FO');
%save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half','.mat'),'FC_ICA_FO','FC_n_FO','FC_ICA_PO','FC_n_PO','FC_ICA_FP','FC_n_FP','FC_SQ_FO','FC_SQ_PO','FC_SQ_FP');

end

function FC_groupICA_ROI_grouped_by_lobe_clusters(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

[FC_ICA_cluster_FO,FC_n_cluster_FO] = getFC_ICA_cluster(FC_ICA_FO,FC_n_FO);
[FC_ICA_cluster_PO,FC_n_cluster_PO] = getFC_ICA_cluster(FC_ICA_PO,FC_n_PO);
[FC_ICA_cluster_FP,FC_n_cluster_FP] = getFC_ICA_cluster(FC_ICA_FP,FC_n_FP);

[FC_SQ_cluster_FO] = getFC_SQ_cluster(FC_ICA_cluster_FO,FC_n_cluster_FO);
[FC_SQ_cluster_PO] = getFC_SQ_cluster(FC_ICA_cluster_PO,FC_n_cluster_PO);
[FC_SQ_cluster_FP] = getFC_SQ_cluster(FC_ICA_cluster_FP,FC_n_cluster_FP);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe-cluster','.mat'),'FC_ICA_cluster_FO','FC_n_cluster_FO','FC_ICA_cluster_PO','FC_n_cluster_PO','FC_ICA_cluster_FP','FC_n_cluster_FP','FC_SQ_cluster_FO','FC_SQ_cluster_PO','FC_SQ_cluster_FP');

end

function [FC_ICA, FC_n] = getFC_ICA(settings,idx_area1,idx_area2,area1_label,area2_label)

matlabpool('OPEN','AttachedFiles','C:\Program Files\MATLAB\R2014a\toolbox\matlab\specfun\betainc.m');
    
ICA_preprocessed = 'groupICA';
nTR = 300;
sample = 10;

area_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area');

FC_n.area1_n = length(idx_area1);
FC_n.area2_n = length(idx_area2);

FC_n.area1_label = area1_label;
FC_n.area2_label = area2_label;

for irun=1:2
    
    all_area1 = [];
    all_area2 = [];
    
    all_run = [];
    
    nICs = 0;
    
    for iROI=1:FC_n.area1_n
   
        area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(irun),'-','aal2std_',int2str(idx_area1(iROI)),'\','melodic_mix'));
           
        nICs = nICs + size(area_ic,2);
        %nICs = nICs + sample;
        
        FC_n.run(irun).nIC_area1(iROI) = size(area_ic,2);
        %FC_n.run(irun).nIC_area1(iROI) = sample;
        
        all_area1 = [all_area1,area_ic];
        %all_area1 = [all_area1,area_ic(:,1:sample)];
        
        clear area_ic;
        
    end
    
    FC_n.run(irun).nIC_area1_total = nICs;
    
    nICs = 0;
    
    for iROI=1:FC_n.area2_n
   
        area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(irun),'-','aal2std_',int2str(idx_area2(iROI)),'\','melodic_mix'));
           
        nICs = nICs + size(area_ic,2);
        %nICs = nICs + sample;
        
        FC_n.run(irun).nIC_area2(iROI) = size(area_ic,2);
        %FC_n.run(irun).nIC_area2(iROI) = sample;
        
        all_area2 = [all_area2,area_ic];
        %all_area2 = [all_area2,area_ic(:,1:sample)];
        
        clear area_ic;
        
    end
    
    FC_n.run(irun).nIC_area2_total = nICs;
    
    all_run = [all_area1,all_area2];
    
    high_run_ic = all_run(1:nTR,:);
    rest_run_ic = all_run(nTR+1:end,:);
    
    n_ic_high_run = size(high_run_ic,2);
    n_ic_rest_run = size(rest_run_ic,2);
    
    disp(strcat('doing high:',area1_label,'-',area2_label));
    disp(strcat('ic:',int2str(n_ic_high_run)));
    
    %[rho_high, pval_high] = corr(high_run_ic);
    
    rho_high = zeros(n_ic_high_run,n_ic_high_run);
    pval_high = zeros(n_ic_high_run,n_ic_high_run);
    
    for iic_line=1:n_ic_high_run
        
        X = high_run_ic(:,iic_line);
        
        %gpu_X = gpuArray(X);
        
        rho_col = [];
        pval_col = [];
        
        columns = iic_line:n_ic_high_run;
       
        %tic
        
        parfor iic_col=1:length(columns)
            
            Y = high_run_ic(:,columns(iic_col));
            
            %gpu_Y = gpuArray(Y);
            
            [rho,pval] = corrcoef(X,Y);
            %rho = corr2(gpu_X,gpu_Y);
            
%             t = (gather(rho)*sqrt(nTR-2))/sqrt(1-gather(rho)^2);
%             pval = 1-tcdf(t,nTR-1);
            
            rho_col(iic_col) = rho(1,2);
            pval_col(iic_col) = pval(1,2);
            
            %rho_col(iic_col) = gather(rho);
            %pval_col(iic_col) = pval;
    
        end
        
        %toc
        
        rho_high(iic_line,iic_line:n_ic_high_run) = rho_col;
        pval_high(iic_line,iic_line:n_ic_high_run) = pval_col;
        
        rho_high(iic_line:n_ic_high_run,iic_line) = rho_col';
        pval_high(iic_line:n_ic_high_run,iic_line) = pval_col';
        
    end
    
    FC_ICA.run(irun).rho_high = rho_high;
    FC_ICA.run(irun).pval_high = pval_high;
    
    disp(strcat('doing rest:',area1_label,'-',area2_label));
    disp(strcat('ic:',int2str(n_ic_rest_run)));
    
    %[rho_rest, pval_rest] = corr(rest_run_ic);
    
    rho_rest = zeros(n_ic_rest_run,n_ic_rest_run);
    pval_rest = zeros(n_ic_rest_run,n_ic_rest_run);
    
    for iic_line=1:n_ic_rest_run
        
        X = rest_run_ic(:,iic_line);
        
        %gpu_X = gpuArray(X);
        
        rho_col = [];
        pval_col = [];
        
        columns = iic_line:n_ic_rest_run;
       
        %tic

        parfor iic_col=1:length(columns)
            
            Y = rest_run_ic(:,columns(iic_col));
             
            %gpu_Y = gpuArray(Y);
            
            [rho,pval] = corrcoef(X,Y);
            %rho = corr2(gpu_X,gpu_Y);
            
            %t = (gather(rho)*sqrt(nTR-2))/sqrt(1-gather(rho)^2);
            %pval = 1-tcdf(t,nTR-1);
            
            rho_col(iic_col) = rho(1,2);
            pval_col(iic_col) = pval(1,2);
            
            %rho_col(iic_col) = gather(rho);
            %pval_col(iic_col) = pval;
    
        end
        
        %toc

        rho_rest(iic_line,iic_line:n_ic_rest_run) = rho_col;
        pval_rest(iic_line,iic_line:n_ic_rest_run) = pval_col;
        
        rho_rest(iic_line:n_ic_rest_run,iic_line) = rho_col';
        pval_rest(iic_line:n_ic_rest_run,iic_line) = pval_col';
        
    end
    
    FC_ICA.run(irun).rho_rest = rho_rest;
    FC_ICA.run(irun).pval_rest = pval_rest;
    
    clear high_run_ic;
    clear rest_run_ic;
    
end

matlabpool CLOSE

end

function [FC_ICA, FC_n] = getFC_ICA_Half(settings,idx_area1,idx_area2,area1_label,area2_label)

matlabpool('OPEN','AttachedFiles','C:\Program Files\MATLAB\R2014a\toolbox\matlab\specfun\betainc.m');
    
ICA_preprocessed = 'groupICA';
nTR = 300;
sample = 10;

area_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area');

FC_n.area1_n = length(idx_area1);
FC_n.area2_n = length(idx_area2);

FC_n.area1_label = area1_label;
FC_n.area2_label = area2_label;

for irun=1:2
    
    all_area1 = [];
    all_area2 = [];
    
    all_run = [];
    
    nICs = 0;
    
    for iROI=1:FC_n.area1_n
   
        area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(irun),'-','aal2std_',int2str(idx_area1(iROI)),'\','melodic_mix'));
           
        nICs = nICs + size(area_ic,2);
        %nICs = nICs + sample;
        
        FC_n.run(irun).nIC_area1(iROI) = size(area_ic,2);
        %FC_n.run(irun).nIC_area1(iROI) = sample;
        
        all_area1 = [all_area1,area_ic];
        %all_area1 = [all_area1,area_ic(:,1:sample)];
        
        clear area_ic;
        
    end
    
    FC_n.run(irun).nIC_area1_total = nICs;
    
    nICs = 0;
    
    for iROI=1:FC_n.area2_n
   
        area_ic = load(strcat(area_preprocessed_folder,'\','high-rest','\','high-rest-run-',int2str(irun),'-','aal2std_',int2str(idx_area2(iROI)),'\','melodic_mix'));
           
        nICs = nICs + size(area_ic,2);
        %nICs = nICs + sample;
        
        FC_n.run(irun).nIC_area2(iROI) = size(area_ic,2);
        %FC_n.run(irun).nIC_area2(iROI) = sample;
        
        all_area2 = [all_area2,area_ic];
        %all_area2 = [all_area2,area_ic(:,1:sample)];
        
        clear area_ic;
        
    end
    
    FC_n.run(irun).nIC_area2_total = nICs;
    
    all_run = [all_area1,all_area2];
    
    high_run_ic = all_run(1:nTR,:);
    rest_run_ic = all_run(nTR+1:end,:);
    
    n_ic_high_run = size(high_run_ic,2);
    n_ic_rest_run = size(rest_run_ic,2);
    
    disp(strcat('doing high:',area1_label,'-',area2_label));
    disp(strcat('ic:',int2str(n_ic_high_run)));
    
    %[rho_high, pval_high] = corr(high_run_ic);
    
    rho_high_FH = zeros(n_ic_high_run,n_ic_high_run);
    pval_high_FH = zeros(n_ic_high_run,n_ic_high_run);
    
    rho_high_SH = zeros(n_ic_high_run,n_ic_high_run);
    pval_high_SH = zeros(n_ic_high_run,n_ic_high_run);
    
    for iic_line=1:n_ic_high_run
        
        X_FH = high_run_ic(1:(nTR/2),iic_line);
        X_SH = high_run_ic(((nTR/2)+1):end,iic_line);
        
        %gpu_X = gpuArray(X);
        
        rho_col_FH = [];
        pval_col_FH = [];
        
        rho_col_SH = [];
        pval_col_SH = [];
        
        columns = iic_line:n_ic_high_run;
       
        %tic
        
        parfor iic_col=1:length(columns)
            
            Y_FH = high_run_ic(1:(nTR/2),columns(iic_col));
            Y_SH = high_run_ic(((nTR/2)+1):end,columns(iic_col));
            
            %gpu_Y = gpuArray(Y);
            
            [rho_FH,pval_FH] = corrcoef(X_FH,Y_FH);
            [rho_SH,pval_SH] = corrcoef(X_SH,Y_SH);
            %rho = corr2(gpu_X,gpu_Y);
            
%             t = (gather(rho)*sqrt(nTR-2))/sqrt(1-gather(rho)^2);
%             pval = 1-tcdf(t,nTR-1);
            
            rho_col_FH(iic_col) = rho_FH(1,2);
            pval_col_FH(iic_col) = pval_FH(1,2);
            
            rho_col_SH(iic_col) = rho_SH(1,2);
            pval_col_SH(iic_col) = pval_SH(1,2);
            
            %rho_col(iic_col) = gather(rho);
            %pval_col(iic_col) = pval;
    
        end
        
        %toc
        
        rho_high_FH(iic_line,iic_line:n_ic_high_run) = rho_col_FH;
        pval_high_FH(iic_line,iic_line:n_ic_high_run) = pval_col_FH;
        
        rho_high_FH(iic_line:n_ic_high_run,iic_line) = rho_col_FH';
        pval_high_FH(iic_line:n_ic_high_run,iic_line) = pval_col_FH';
        
        rho_high_SH(iic_line,iic_line:n_ic_high_run) = rho_col_SH;
        pval_high_SH(iic_line,iic_line:n_ic_high_run) = pval_col_SH;
        
        rho_high_SH(iic_line:n_ic_high_run,iic_line) = rho_col_SH';
        pval_high_SH(iic_line:n_ic_high_run,iic_line) = pval_col_SH';
        
    end
    
    FC_ICA.run(irun).rho_high_FH = rho_high_FH;
    FC_ICA.run(irun).pval_high_FH = pval_high_FH;
    
    FC_ICA.run(irun).rho_high_SH = rho_high_SH;
    FC_ICA.run(irun).pval_high_SH = pval_high_SH;
    
    disp(strcat('doing rest:',area1_label,'-',area2_label));
    disp(strcat('ic:',int2str(n_ic_rest_run)));
    
    %[rho_rest, pval_rest] = corr(rest_run_ic);
    
    rho_rest_FH = zeros(n_ic_rest_run,n_ic_rest_run);
    pval_rest_FH = zeros(n_ic_rest_run,n_ic_rest_run);
    
    rho_rest_SH = zeros(n_ic_rest_run,n_ic_rest_run);
    pval_rest_SH = zeros(n_ic_rest_run,n_ic_rest_run);
    
    for iic_line=1:n_ic_rest_run
        
        X_FH = rest_run_ic(1:(nTR/2),iic_line);
        X_SH = rest_run_ic(((nTR/2)+1):end,iic_line);
        
        %gpu_X = gpuArray(X);
        
        rho_col_FH = [];
        pval_col_FH = [];
        
        rho_col_SH = [];
        pval_col_SH = [];
        
        columns = iic_line:n_ic_rest_run;
       
        %tic

        parfor iic_col=1:length(columns)
            
            Y_FH = rest_run_ic(1:(nTR/2),columns(iic_col));
            Y_SH = rest_run_ic(((nTR/2)+1):end,columns(iic_col));
             
            %gpu_Y = gpuArray(Y);
            
            [rho_FH,pval_FH] = corrcoef(X_FH,Y_FH);
            [rho_SH,pval_SH] = corrcoef(X_SH,Y_SH);
            %rho = corr2(gpu_X,gpu_Y);
            
            %t = (gather(rho)*sqrt(nTR-2))/sqrt(1-gather(rho)^2);
            %pval = 1-tcdf(t,nTR-1);
            
            rho_col_FH(iic_col) = rho_FH(1,2);
            pval_col_FH(iic_col) = pval_FH(1,2);
            
            rho_col_SH(iic_col) = rho_SH(1,2);
            pval_col_SH(iic_col) = pval_SH(1,2);
            
            %rho_col(iic_col) = gather(rho);
            %pval_col(iic_col) = pval;
    
        end
        
        %toc

        rho_rest_FH(iic_line,iic_line:n_ic_rest_run) = rho_col_FH;
        pval_rest_FH(iic_line,iic_line:n_ic_rest_run) = pval_col_FH;
        
        rho_rest_FH(iic_line:n_ic_rest_run,iic_line) = rho_col_FH';
        pval_rest_FH(iic_line:n_ic_rest_run,iic_line) = pval_col_FH';
        
        rho_rest_SH(iic_line,iic_line:n_ic_rest_run) = rho_col_SH;
        pval_rest_SH(iic_line,iic_line:n_ic_rest_run) = pval_col_SH;
        
        rho_rest_SH(iic_line:n_ic_rest_run,iic_line) = rho_col_SH';
        pval_rest_SH(iic_line:n_ic_rest_run,iic_line) = pval_col_SH';
        
    end
    
    FC_ICA.run(irun).rho_rest_FH = rho_rest_FH;
    FC_ICA.run(irun).pval_rest_FH = pval_rest_FH;
    
    FC_ICA.run(irun).rho_rest_SH = rho_rest_SH;
    FC_ICA.run(irun).pval_rest_SH = pval_rest_SH;
    
    clear high_run_ic;
    clear rest_run_ic;
    
end

matlabpool CLOSE

end

function [FC_ICA_cluster,FC_n_cluster] = getFC_ICA_cluster(FC_ICA,FC_n)

num_Clust = 8;
[T_km,C,sumd,D] = kmeans(CorMx, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);


end

function [FC_SQ] = getFC_SQ(FC_ICA,FC_n)

for irun=1:2

    rho_high_ic = FC_ICA.run(irun).rho_high;
    rho_rest_ic = FC_ICA.run(irun).rho_rest;
    
    rho_high_roi = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    rho_rest_roi = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    
    for iroi=1:FC_n.area1_n
       
        if iroi==1
            nIC_start_line = 1;
        else
            nIC_start_line = sum(FC_n.run(irun).nIC_area1(1:iroi-1)) + 1;  
        end
                
        nIC_end_line = sum(FC_n.run(irun).nIC_area1(1:iroi));
            
        for iiroi1=1:FC_n.area1_n
           
            if iiroi1==1
                nIC_start = 1;
            else
                nIC_start = sum(FC_n.run(irun).nIC_area1(1:iiroi1-1)) + 1;
            end
                
            nIC_end = sum(FC_n.run(irun).nIC_area1(1:iiroi1));
            
            rho_high_roi(iroi,iiroi1) = sumsqr(rho_high_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi(iroi,iiroi1) = sumsqr(rho_rest_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
        for iiroi2=1:FC_n.area2_n
           
            if iiroi2==1
                nIC_start = FC_n.run(irun).nIC_area1_total + 1;
            else
                nIC_start = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2-1)) + 1;
            end
                
            nIC_end = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2));
            
            rho_high_roi(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
    end
    
    for iroi=1:FC_n.area2_n
       
        if iroi==1
            nIC_start_line = FC_n.run(irun).nIC_area1_total + 1;
        else
            nIC_start_line = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iroi-1)) + 1;  
        end
                
        nIC_end_line = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iroi));
            
        for iiroi1=1:FC_n.area1_n
           
            if iiroi1==1
                nIC_start = 1;
            else
                nIC_start = sum(FC_n.run(irun).nIC_area1(1:iiroi1-1)) + 1;
            end
                
            nIC_end = sum(FC_n.run(irun).nIC_area1(1:iiroi1));
            
            rho_high_roi(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_high_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_rest_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
        for iiroi2=1:FC_n.area2_n
           
            if iiroi2==1
                nIC_start = FC_n.run(irun).nIC_area1_total + 1;
            else
                nIC_start = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2-1)) + 1;
            end
                
            nIC_end = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2));
            
            rho_high_roi(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
       end
       
       FC_SQ.run(irun).rho_high_roi = rho_high_roi;
       FC_SQ.run(irun).rho_rest_roi = rho_rest_roi;
        
end
   

end

function [FC_SQ] = getFC_SQ_Half(FC_ICA,FC_n)

for irun=1:2

    rho_high_ic_FH = FC_ICA.run(irun).rho_high_FH;
    rho_rest_ic_FH = FC_ICA.run(irun).rho_rest_FH;
    
    rho_high_ic_SH = FC_ICA.run(irun).rho_high_SH;
    rho_rest_ic_SH = FC_ICA.run(irun).rho_rest_SH;
    
    rho_high_roi_FH = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    rho_rest_roi_FH = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    
    rho_high_roi_SH = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    rho_rest_roi_SH = zeros(FC_n.area1_n+FC_n.area2_n,FC_n.area1_n+FC_n.area2_n);
    
    for iroi=1:FC_n.area1_n
       
        if iroi==1
            nIC_start_line = 1;
        else
            nIC_start_line = sum(FC_n.run(irun).nIC_area1(1:iroi-1)) + 1;  
        end
                
        nIC_end_line = sum(FC_n.run(irun).nIC_area1(1:iroi));
            
        for iiroi1=1:FC_n.area1_n
           
            if iiroi1==1
                nIC_start = 1;
            else
                nIC_start = sum(FC_n.run(irun).nIC_area1(1:iiroi1-1)) + 1;
            end
                
            nIC_end = sum(FC_n.run(irun).nIC_area1(1:iiroi1));
            
            rho_high_roi_FH(iroi,iiroi1) = sumsqr(rho_high_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_FH(iroi,iiroi1) = sumsqr(rho_rest_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
            rho_high_roi_SH(iroi,iiroi1) = sumsqr(rho_high_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_SH(iroi,iiroi1) = sumsqr(rho_rest_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
        for iiroi2=1:FC_n.area2_n
           
            if iiroi2==1
                nIC_start = FC_n.run(irun).nIC_area1_total + 1;
            else
                nIC_start = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2-1)) + 1;
            end
                
            nIC_end = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2));
            
            rho_high_roi_FH(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_FH(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
            rho_high_roi_SH(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_SH(iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
    end
    
    for iroi=1:FC_n.area2_n
       
        if iroi==1
            nIC_start_line = FC_n.run(irun).nIC_area1_total + 1;
        else
            nIC_start_line = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iroi-1)) + 1;  
        end
                
        nIC_end_line = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iroi));
            
        for iiroi1=1:FC_n.area1_n
           
            if iiroi1==1
                nIC_start = 1;
            else
                nIC_start = sum(FC_n.run(irun).nIC_area1(1:iiroi1-1)) + 1;
            end
                
            nIC_end = sum(FC_n.run(irun).nIC_area1(1:iiroi1));
            
            rho_high_roi_FH(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_high_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_FH(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_rest_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
            rho_high_roi_SH(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_high_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_SH(FC_n.area1_n+iroi,iiroi1) = sumsqr(rho_rest_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
        for iiroi2=1:FC_n.area2_n
           
            if iiroi2==1
                nIC_start = FC_n.run(irun).nIC_area1_total + 1;
            else
                nIC_start = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2-1)) + 1;
            end
                
            nIC_end = FC_n.run(irun).nIC_area1_total + sum(FC_n.run(irun).nIC_area2(1:iiroi2));
            
            rho_high_roi_FH(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_FH(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic_FH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
            rho_high_roi_SH(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_high_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            rho_rest_roi_SH(FC_n.area1_n+iroi,FC_n.area1_n+iiroi2) = sumsqr(rho_rest_ic_SH(nIC_start_line:nIC_end_line,nIC_start:nIC_end));
            
        end
        
       end
       
       FC_SQ.run(irun).rho_high_roi_FH = rho_high_roi_FH;
       FC_SQ.run(irun).rho_rest_roi_FH = rho_rest_roi_FH;
       
       FC_SQ.run(irun).rho_high_roi_SH = rho_high_roi_SH;
       FC_SQ.run(irun).rho_rest_roi_SH = rho_rest_roi_SH;
        
end
   
end

function [FC_SQ_cluster] = getFC_SQ_cluster(FC_ICA_cluster,FC_n_cluster)



end

function plot_FC_ICs(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

ViewFC(settings,1,FC_n_FO.area1_label, FC_n_FO.area2_label, FC_ICA_FO.run(1).rho_high, FC_ICA_FO.run(1).pval_high, FC_ICA_FO.run(1).rho_rest, FC_ICA_FO.run(1).pval_rest, FC_n_FO.run(1).nIC_area1_total, FC_n_FO.run(1).nIC_area2_total );
ViewFC(settings,2,FC_n_FO.area1_label, FC_n_FO.area2_label, FC_ICA_FO.run(2).rho_high, FC_ICA_FO.run(2).pval_high, FC_ICA_FO.run(2).rho_rest, FC_ICA_FO.run(2).pval_rest, FC_n_FO.run(2).nIC_area1_total, FC_n_FO.run(2).nIC_area2_total );

ViewFC(settings,1,FC_n_FP.area1_label, FC_n_FP.area2_label, FC_ICA_FP.run(1).rho_high, FC_ICA_FP.run(1).pval_high, FC_ICA_FP.run(1).rho_rest, FC_ICA_FP.run(1).pval_rest, FC_n_FP.run(1).nIC_area1_total, FC_n_FP.run(1).nIC_area2_total );
ViewFC(settings,2,FC_n_FP.area1_label, FC_n_FP.area2_label, FC_ICA_FP.run(2).rho_high, FC_ICA_FP.run(2).pval_high, FC_ICA_FP.run(2).rho_rest, FC_ICA_FP.run(2).pval_rest, FC_n_FP.run(2).nIC_area1_total, FC_n_FP.run(2).nIC_area2_total );

ViewFC(settings,1,FC_n_PO.area1_label, FC_n_PO.area2_label, FC_ICA_PO.run(1).rho_high, FC_ICA_PO.run(1).pval_high, FC_ICA_PO.run(1).rho_rest, FC_ICA_PO.run(1).pval_rest, FC_n_PO.run(1).nIC_area1_total, FC_n_PO.run(1).nIC_area2_total );
ViewFC(settings,2,FC_n_PO.area1_label, FC_n_PO.area2_label, FC_ICA_PO.run(2).rho_high, FC_ICA_PO.run(2).pval_high, FC_ICA_PO.run(2).rho_rest, FC_ICA_PO.run(2).pval_rest, FC_n_PO.run(2).nIC_area1_total, FC_n_PO.run(2).nIC_area2_total );

end

function plot_FC_SumOfSquare(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

%%% HIGH
rho_eye = eye(size(FC_SQ_FO.run(1).rho_high_roi));
FC_SQ_FO.run(1).rho_high_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FO.run(2).rho_high_roi));
FC_SQ_FO.run(2).rho_high_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FP.run(1).rho_high_roi));
FC_SQ_FP.run(1).rho_high_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FP.run(2).rho_high_roi));
FC_SQ_FP.run(2).rho_high_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_PO.run(1).rho_high_roi));
FC_SQ_PO.run(1).rho_high_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_PO.run(2).rho_high_roi));
FC_SQ_PO.run(2).rho_high_roi(logical(rho_eye)) = 0;

%% REST

rho_eye = eye(size(FC_SQ_FO.run(1).rho_rest_roi));
FC_SQ_FO.run(1).rho_rest_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FO.run(2).rho_rest_roi));
FC_SQ_FO.run(2).rho_rest_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FP.run(1).rho_rest_roi));
FC_SQ_FP.run(1).rho_rest_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_FP.run(2).rho_rest_roi));
FC_SQ_FP.run(2).rho_rest_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_PO.run(1).rho_rest_roi));
FC_SQ_PO.run(1).rho_rest_roi(logical(rho_eye)) = 0;

rho_eye = eye(size(FC_SQ_PO.run(2).rho_rest_roi));
FC_SQ_PO.run(2).rho_rest_roi(logical(rho_eye)) = 0;

max_high = max([max(FC_SQ_FO.run(1).rho_high_roi(:)), max(FC_SQ_FO.run(2).rho_high_roi(:)), max(FC_SQ_FP.run(1).rho_high_roi(:)), max(FC_SQ_FP.run(2).rho_high_roi(:)), max(FC_SQ_PO.run(1).rho_high_roi(:)), max(FC_SQ_PO.run(2).rho_high_roi(:))]);
min_high = min([min(FC_SQ_FO.run(1).rho_high_roi(:)), min(FC_SQ_FO.run(2).rho_high_roi(:)), min(FC_SQ_FP.run(1).rho_high_roi(:)), min(FC_SQ_FP.run(2).rho_high_roi(:)), min(FC_SQ_PO.run(1).rho_high_roi(:)), min(FC_SQ_PO.run(2).rho_high_roi(:))]);

max_rest = max([max(FC_SQ_FO.run(1).rho_rest_roi(:)), max(FC_SQ_FO.run(2).rho_rest_roi(:)), max(FC_SQ_FP.run(1).rho_rest_roi(:)), max(FC_SQ_FP.run(2).rho_rest_roi(:)), max(FC_SQ_PO.run(1).rho_rest_roi(:)), max(FC_SQ_PO.run(2).rho_rest_roi(:))]);
min_rest = min([min(FC_SQ_FO.run(1).rho_rest_roi(:)), min(FC_SQ_FO.run(2).rho_rest_roi(:)), min(FC_SQ_FP.run(1).rho_rest_roi(:)), min(FC_SQ_FP.run(2).rho_rest_roi(:)), min(FC_SQ_PO.run(1).rho_rest_roi(:)), min(FC_SQ_PO.run(2).rho_rest_roi(:))]);

max_all = max(max_high,max_rest);
min_all = min(min_high,min_rest);

CLIM = [ceil(min_all) ceil(max_all)];
CLIM = [0 1000];

ViewFC_SQ(settings,1,FC_n_FO.area1_label, FC_n_FO.area2_label, FC_SQ_FO.run(1).rho_high_roi, FC_SQ_FO.run(1).rho_rest_roi, FC_n_FO.area1_n, FC_n_FO.area2_n, CLIM );
ViewFC_SQ(settings,2,FC_n_FO.area1_label, FC_n_FO.area2_label, FC_SQ_FO.run(2).rho_high_roi, FC_SQ_FO.run(2).rho_rest_roi, FC_n_FO.area1_n, FC_n_FO.area2_n, CLIM );

ViewFC_SQ(settings,1,FC_n_FP.area1_label, FC_n_FP.area2_label, FC_SQ_FP.run(1).rho_high_roi, FC_SQ_FP.run(1).rho_rest_roi, FC_n_FP.area1_n, FC_n_FP.area2_n, CLIM );
ViewFC_SQ(settings,2,FC_n_FP.area1_label, FC_n_FP.area2_label, FC_SQ_FP.run(2).rho_high_roi, FC_SQ_FP.run(2).rho_rest_roi, FC_n_FP.area1_n, FC_n_FP.area2_n, CLIM );

ViewFC_SQ(settings,1,FC_n_PO.area1_label, FC_n_PO.area2_label, FC_SQ_PO.run(1).rho_high_roi, FC_SQ_PO.run(1).rho_rest_roi, FC_n_PO.area1_n, FC_n_PO.area2_n, CLIM );
ViewFC_SQ(settings,2,FC_n_PO.area1_label, FC_n_PO.area2_label, FC_SQ_PO.run(2).rho_high_roi, FC_SQ_PO.run(2).rho_rest_roi, FC_n_PO.area1_n, FC_n_PO.area2_n, CLIM );


end

function plot_FC_SumOfSquare_Together(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

F_n = FC_n_FO.area1_n;
O_n = FC_n_FO.area2_n;
P_n = FC_n_PO.area1_n;

%%% HIGH ATTENTION

FF_run1 = FC_SQ_FO.run(1).rho_high_roi(1:F_n,1:F_n);
FF_run2 = FC_SQ_FO.run(2).rho_high_roi(1:F_n,1:F_n);

OO_run1 = FC_SQ_FO.run(1).rho_high_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
OO_run2 = FC_SQ_FO.run(2).rho_high_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);

PP_run1 = FC_SQ_PO.run(1).rho_high_roi(1:P_n,1:P_n);
PP_run2 = FC_SQ_PO.run(2).rho_high_roi(1:P_n,1:P_n);

FO_run1 = FC_SQ_FO.run(1).rho_high_roi(1:F_n,F_n+1:F_n+O_n);
FO_run2 = FC_SQ_FO.run(2).rho_high_roi(1:F_n,F_n+1:F_n+O_n);

PO_run1 = FC_SQ_PO.run(1).rho_high_roi(1:P_n,P_n+1:P_n+O_n);
PO_run2 = FC_SQ_PO.run(2).rho_high_roi(1:P_n,P_n+1:P_n+O_n);

FP_run1 = FC_SQ_FP.run(1).rho_high_roi(1:F_n,F_n+1:F_n+P_n);
FP_run2 = FC_SQ_FP.run(2).rho_high_roi(1:F_n,F_n+1:F_n+P_n);

F_run1 = [FF_run1,FP_run1,FO_run1];
F_run2 = [FF_run2,FP_run2,FO_run2];

P_run1 = [FP_run1',PP_run1,PO_run1];
P_run2 = [FP_run2',PP_run2,PO_run2];

O_run1 = [FO_run1',PO_run1',OO_run1];
O_run2 = [FO_run2',PO_run2',OO_run2];

AAL_run1 = [F_run1;P_run1;O_run1];
AAL_run2 = [F_run2;P_run2;O_run2];

correlation = corr([AAL_run1(:),AAL_run2(:)]);
disp(strcat('Correlation btw. Runs - High Attention:',num2str(correlation)));

f = figure;

min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 
 
%max_aal = mean(mean(AAL_run1)) + 2*std(mean(AAL_run1));
max_aal = 300;

CLIM = [min_aal max_aal];

imagesc(AAL_run1,CLIM);

title('High Attention - Run 1');

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

hold on;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1','.pdf'));

f = figure;

min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 
 
%CLIM = [min_aal max_aal];

imagesc(AAL_run2,CLIM);

title('High Attention - Run 2');

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

hold on;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2','.pdf'));

%%% RESTING STATE

FF_run1 = FC_SQ_FO.run(1).rho_rest_roi(1:F_n,1:F_n);
FF_run2 = FC_SQ_FO.run(2).rho_rest_roi(1:F_n,1:F_n);

OO_run1 = FC_SQ_FO.run(1).rho_rest_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
OO_run2 = FC_SQ_FO.run(2).rho_rest_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);

PP_run1 = FC_SQ_PO.run(1).rho_rest_roi(1:P_n,1:P_n);
PP_run2 = FC_SQ_PO.run(2).rho_rest_roi(1:P_n,1:P_n);

FO_run1 = FC_SQ_FO.run(1).rho_rest_roi(1:F_n,F_n+1:F_n+O_n);
FO_run2 = FC_SQ_FO.run(2).rho_rest_roi(1:F_n,F_n+1:F_n+O_n);

PO_run1 = FC_SQ_PO.run(1).rho_rest_roi(1:P_n,P_n+1:P_n+O_n);
PO_run2 = FC_SQ_PO.run(2).rho_rest_roi(1:P_n,P_n+1:P_n+O_n);

FP_run1 = FC_SQ_FP.run(1).rho_rest_roi(1:F_n,F_n+1:F_n+P_n);
FP_run2 = FC_SQ_FP.run(2).rho_rest_roi(1:F_n,F_n+1:F_n+P_n);

F_run1 = [FF_run1,FP_run1,FO_run1];
F_run2 = [FF_run2,FP_run2,FO_run2];

P_run1 = [FP_run1',PP_run1,PO_run1];
P_run2 = [FP_run2',PP_run2,PO_run2];

O_run1 = [FO_run1',PO_run1',OO_run1];
O_run2 = [FO_run2',PO_run2',OO_run2];

AAL_run1 = [F_run1;P_run1;O_run1];
AAL_run2 = [F_run2;P_run2;O_run2];

correlation = corr([AAL_run1(:),AAL_run2(:)]);
disp(strcat('Correlation btw. Runs - Resting State:',num2str(correlation)));

f = figure;

min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 
 
%CLIM = [min_aal max_aal];

imagesc(AAL_run1,CLIM);

title('Resting State - Run 1');

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

hold on;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1','.pdf'));

f = figure;

min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 
 
%CLIM = [min_aal max_aal];

imagesc(AAL_run2,CLIM);

title('Resting State - Run 2');

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

hold on;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2','.pdf'));


end

function plot_FC_SumOfSquare_Together_Half(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FO','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FP','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-PO','.mat'));


for iFPO=1:2
    
    if iFPO==1
        
        which_one = 'FH';
        
        FC_SQ_FO.run(1).rho_high_roi = FC_SQ_FO.run(1).rho_high_roi_FH;
        FC_SQ_FO.run(2).rho_high_roi = FC_SQ_FO.run(2).rho_high_roi_FH;
    
        FC_SQ_FO.run(1).rho_rest_roi = FC_SQ_FO.run(1).rho_rest_roi_FH;
        FC_SQ_FO.run(2).rho_rest_roi = FC_SQ_FO.run(2).rho_rest_roi_FH;
        
        FC_SQ_PO.run(1).rho_high_roi = FC_SQ_PO.run(1).rho_high_roi_FH;
        FC_SQ_PO.run(2).rho_high_roi = FC_SQ_PO.run(2).rho_high_roi_FH;
    
        FC_SQ_PO.run(1).rho_rest_roi = FC_SQ_PO.run(1).rho_rest_roi_FH;
        FC_SQ_PO.run(2).rho_rest_roi = FC_SQ_PO.run(2).rho_rest_roi_FH;
        
        FC_SQ_FP.run(1).rho_high_roi = FC_SQ_FP.run(1).rho_high_roi_FH;
        FC_SQ_FP.run(2).rho_high_roi = FC_SQ_FP.run(2).rho_high_roi_FH;
    
        FC_SQ_FP.run(1).rho_rest_roi = FC_SQ_FP.run(1).rho_rest_roi_FH;
        FC_SQ_FP.run(2).rho_rest_roi = FC_SQ_FP.run(2).rho_rest_roi_FH;
        
    else
        
        which_one = 'SH';
        
        FC_SQ_FO.run(1).rho_high_roi = FC_SQ_FO.run(1).rho_high_roi_SH;
        FC_SQ_FO.run(2).rho_high_roi = FC_SQ_FO.run(2).rho_high_roi_SH;
    
        FC_SQ_FO.run(1).rho_rest_roi = FC_SQ_FO.run(1).rho_rest_roi_SH;
        FC_SQ_FO.run(2).rho_rest_roi = FC_SQ_FO.run(2).rho_rest_roi_SH;
        
        FC_SQ_PO.run(1).rho_high_roi = FC_SQ_PO.run(1).rho_high_roi_SH;
        FC_SQ_PO.run(2).rho_high_roi = FC_SQ_PO.run(2).rho_high_roi_SH;
    
        FC_SQ_PO.run(1).rho_rest_roi = FC_SQ_PO.run(1).rho_rest_roi_SH;
        FC_SQ_PO.run(2).rho_rest_roi = FC_SQ_PO.run(2).rho_rest_roi_SH;
        
        FC_SQ_FP.run(1).rho_high_roi = FC_SQ_FP.run(1).rho_high_roi_SH;
        FC_SQ_FP.run(2).rho_high_roi = FC_SQ_FP.run(2).rho_high_roi_SH;
    
        FC_SQ_FP.run(1).rho_rest_roi = FC_SQ_FP.run(1).rho_rest_roi_SH;
        FC_SQ_FP.run(2).rho_rest_roi = FC_SQ_FP.run(2).rho_rest_roi_SH;
        
    end
    
    F_n = FC_n_FO.area1_n;
    O_n = FC_n_FO.area2_n;
    P_n = FC_n_PO.area1_n;

    %%% HIGH ATTENTION

    FF_run1 = FC_SQ_FO.run(1).rho_high_roi(1:F_n,1:F_n);
    FF_run2 = FC_SQ_FO.run(2).rho_high_roi(1:F_n,1:F_n);

    OO_run1 = FC_SQ_FO.run(1).rho_high_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
    OO_run2 = FC_SQ_FO.run(2).rho_high_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);

    PP_run1 = FC_SQ_PO.run(1).rho_high_roi(1:P_n,1:P_n);
    PP_run2 = FC_SQ_PO.run(2).rho_high_roi(1:P_n,1:P_n);

    FO_run1 = FC_SQ_FO.run(1).rho_high_roi(1:F_n,F_n+1:F_n+O_n);
    FO_run2 = FC_SQ_FO.run(2).rho_high_roi(1:F_n,F_n+1:F_n+O_n);

    PO_run1 = FC_SQ_PO.run(1).rho_high_roi(1:P_n,P_n+1:P_n+O_n);
    PO_run2 = FC_SQ_PO.run(2).rho_high_roi(1:P_n,P_n+1:P_n+O_n);

    FP_run1 = FC_SQ_FP.run(1).rho_high_roi(1:F_n,F_n+1:F_n+P_n);
    FP_run2 = FC_SQ_FP.run(2).rho_high_roi(1:F_n,F_n+1:F_n+P_n);

    F_run1 = [FF_run1,FP_run1,FO_run1];
    F_run2 = [FF_run2,FP_run2,FO_run2];

    P_run1 = [FP_run1',PP_run1,PO_run1];
    P_run2 = [FP_run2',PP_run2,PO_run2];

    O_run1 = [FO_run1',PO_run1',OO_run1];
    O_run2 = [FO_run2',PO_run2',OO_run2];

    AAL_run1 = [F_run1;P_run1;O_run1];
    AAL_run2 = [F_run2;P_run2;O_run2];

    correlation = corr([AAL_run1(:),AAL_run2(:)]);
    disp(strcat('Correlation btw. Runs - High Attention:',num2str(correlation)));

    f = figure;

    min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
    max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 

    %max_aal = mean(mean(AAL_run1)) + 2*std(mean(AAL_run1));
    max_aal = 300;

    CLIM = [min_aal max_aal];

    imagesc(AAL_run1,CLIM);

    title(strcat('High Attention - Run 1','-',which_one));

    Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
    TickLabel = { 'F' ; 'P' ; 'O' };
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    hold on;

    line_position = F_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    line_position = F_n + P_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1-',which_one,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1-',which_one,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run1-',which_one,'.pdf'));

    f = figure;

    min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
    max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 

    %CLIM = [min_aal max_aal];

    imagesc(AAL_run2,CLIM);

    title(strcat('High Attention - Run 2','-',which_one));

    Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
    TickLabel = { 'F' ; 'P' ; 'O' };
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    hold on;

    line_position = F_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    line_position = F_n + P_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2-',which_one,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2-',which_one,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','HighAttention-Run2-',which_one,'.pdf'));

    %%% RESTING STATE

    FF_run1 = FC_SQ_FO.run(1).rho_rest_roi(1:F_n,1:F_n);
    FF_run2 = FC_SQ_FO.run(2).rho_rest_roi(1:F_n,1:F_n);

    OO_run1 = FC_SQ_FO.run(1).rho_rest_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
    OO_run2 = FC_SQ_FO.run(2).rho_rest_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);

    PP_run1 = FC_SQ_PO.run(1).rho_rest_roi(1:P_n,1:P_n);
    PP_run2 = FC_SQ_PO.run(2).rho_rest_roi(1:P_n,1:P_n);

    FO_run1 = FC_SQ_FO.run(1).rho_rest_roi(1:F_n,F_n+1:F_n+O_n);
    FO_run2 = FC_SQ_FO.run(2).rho_rest_roi(1:F_n,F_n+1:F_n+O_n);

    PO_run1 = FC_SQ_PO.run(1).rho_rest_roi(1:P_n,P_n+1:P_n+O_n);
    PO_run2 = FC_SQ_PO.run(2).rho_rest_roi(1:P_n,P_n+1:P_n+O_n);

    FP_run1 = FC_SQ_FP.run(1).rho_rest_roi(1:F_n,F_n+1:F_n+P_n);
    FP_run2 = FC_SQ_FP.run(2).rho_rest_roi(1:F_n,F_n+1:F_n+P_n);

    F_run1 = [FF_run1,FP_run1,FO_run1];
    F_run2 = [FF_run2,FP_run2,FO_run2];

    P_run1 = [FP_run1',PP_run1,PO_run1];
    P_run2 = [FP_run2',PP_run2,PO_run2];

    O_run1 = [FO_run1',PO_run1',OO_run1];
    O_run2 = [FO_run2',PO_run2',OO_run2];

    AAL_run1 = [F_run1;P_run1;O_run1];
    AAL_run2 = [F_run2;P_run2;O_run2];

    correlation = corr([AAL_run1(:),AAL_run2(:)]);
    disp(strcat('Correlation btw. Runs - Resting State:',num2str(correlation)));

    f = figure;

    min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
    max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 

    %CLIM = [min_aal max_aal];

    imagesc(AAL_run1,CLIM);

    title(strcat('Resting State - Run 1','-',which_one));

    Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
    TickLabel = { 'F' ; 'P' ; 'O' };
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    hold on;

    line_position = F_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    line_position = F_n + P_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1-',which_one,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1-',which_one,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run1-',which_one,'.pdf'));

    f = figure;

    min_aal = min([min(min(AAL_run1)) min(min(AAL_run2))]); 
    max_aal = max([max(max(AAL_run1)) max(max(AAL_run2))]); 

    %CLIM = [min_aal max_aal];

    imagesc(AAL_run2,CLIM);

    title(strcat('Resting State - Run 2','-',which_one));

    Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
    TickLabel = { 'F' ; 'P' ; 'O' };
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    hold on;

    line_position = F_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    line_position = F_n + P_n;
    plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
    plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2-',which_one,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2-',which_one,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-IC-SQ-AAL-FPO','-','RestingState-Run2-',which_one,'.pdf'));

end
    
end


function ViewFC(settings,run,area1_label, area2_label, rho_att, pval_att, rho_rest, pval_rest, Ncomponent_area1, Ncomponent_area2 )

ROIlabel = strcat(area1_label,'-',area2_label,'-','run-',int2str(run));

fs = 14;

pcriterion = 0.001;
rholimit   = 0.5;

N2 = (length(rho_att(:)))^2;

% prepare correlation plots

rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;

CLIM = [rhomin rhomax];

rhorange = linspace( rhomin, rhomax, 64 );

kk = find( pval_att > pcriterion );
rho_att(kk) = zeros(size(kk));

kk = find( pval_rest > pcriterion );
rho_rest(kk) = zeros(size(kk));


clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

f = figure;
subplot(1,2,1);
hold on;
caxis([rhomin rhomax]);
%h = pcolor( rho_att );
h = imagesc( rho_att, CLIM );
%set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent_area1*[1 1], [0 (Ncomponent_area1+Ncomponent_area2)], 'k', 'LineWidth', 2 );
plot( [0 (Ncomponent_area1+Ncomponent_area2)], 1+ Ncomponent_area1*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent_area1/2 (Ncomponent_area1+Ncomponent_area2/2)];
TickLabel = { area1_label ; area2_label };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('attention', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);


subplot(1,2,2);
hold on;
caxis([rhomin rhomax]);
%h = pcolor( rho_rest );
h = imagesc( rho_rest, CLIM );
%set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent_area1*[1 1], [0 (Ncomponent_area1+Ncomponent_area2)], 'k', 'LineWidth', 2 );
plot( [0 (Ncomponent_area1+Ncomponent_area2)], 1+ Ncomponent_area1*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent_area1/2 (Ncomponent_area1+Ncomponent_area2/2)];
TickLabel = { area1_label ; area2_label };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('resting', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);

suptitle( ROIlabel );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'.pdf'));


end

function ViewFC_SQ(settings,run,area1_label, area2_label, rho_att, rho_rest, Ncomponent_area1, Ncomponent_area2, CLIM )

ROIlabel = strcat(area1_label,'-',area2_label,'-','run-',int2str(run),'-SQ');

rho_att_eye = eye(size(rho_att));
rho_rest_eye = eye(size(rho_rest));

rho_att(logical(rho_att_eye)) = 0;
rho_rest(logical(rho_rest_eye)) = 0;

%CLIM = [0 2];

fs = 14;

clrmp = colormap('jet');
%clrmp(32,:) = [1 1 1];

f = figure;
subplot(1,2,1);
hold on;
caxis(CLIM);
%h = pcolor( rho_att );
h = imagesc( rho_att, CLIM );
%set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( Ncomponent_area1*[1 1], [1 (Ncomponent_area1+Ncomponent_area2)], 'k', 'LineWidth', 2 );
plot( [1 (Ncomponent_area1+Ncomponent_area2)], Ncomponent_area1*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent_area1/2 (Ncomponent_area1+Ncomponent_area2/2)];
TickLabel = { area1_label ; area2_label };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('attention', 'FontSize', fs );


h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);


subplot(1,2,2);
hold on;
caxis(CLIM);
%h = pcolor( rho_rest );
h = imagesc( rho_rest, CLIM );
%set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( Ncomponent_area1*[1 1], [1 (Ncomponent_area1+Ncomponent_area2)], 'k', 'LineWidth', 2 );
plot( [1 (Ncomponent_area1+Ncomponent_area2)], Ncomponent_area1*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent_area1/2 (Ncomponent_area1+Ncomponent_area2/2)];
TickLabel = { area1_label ; area2_label };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('resting', 'FontSize', fs );


h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);

suptitle( ROIlabel );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'-SQ','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'-SQ','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','-',ROIlabel,'-SQ','.pdf'));


end

function plot_FC_lobes(all_settings)

for iset=1:length(all_settings)
    
    settings = all_settings(iset).settings;
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

    F_n = FC_n_FO.area1_n;
    O_n = FC_n_FO.area2_n;
    P_n = FC_n_PO.area1_n;

    for irun=1:2
        
        O(iset).run(irun).rho_high = FC_SQ_FO.run(irun).rho_high_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
        O(iset).run(irun).rho_rest = FC_SQ_FO.run(irun).rho_rest_roi(F_n+1:F_n+O_n,F_n+1:F_n+O_n);
    
        P(iset).run(irun).rho_high = FC_SQ_FP.run(irun).rho_high_roi(F_n+1:F_n+P_n,F_n+1:F_n+P_n);
        P(iset).run(irun).rho_rest = FC_SQ_FP.run(irun).rho_rest_roi(F_n+1:F_n+P_n,F_n+1:F_n+P_n);
    
        F(iset).run(irun).rho_high = FC_SQ_FO.run(irun).rho_high_roi(1:F_n,1:F_n);
        F(iset).run(irun).rho_rest = FC_SQ_FO.run(irun).rho_rest_roi(1:F_n,1:F_n);
    
    end    
  
end

runs = 2;

mean_o_high = 0;
std_o_high = 0;

mean_o_rest = 0;
std_o_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_o_high = mean_o_high + mean(O(iset).run(irun).rho_high(:));
        std_o_high = std_o_high + std(O(iset).run(irun).rho_high(:));
    
        mean_o_rest = mean_o_rest + mean(O(iset).run(irun).rho_rest(:));
        std_o_rest = std_o_rest + std(O(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_o_high = mean_o_high / ( length(all_settings) * runs );
std_o_high = std_o_high / ( length(all_settings) * runs );

mean_o_rest = mean_o_rest / ( length(all_settings) * runs );
std_o_rest = std_o_rest / ( length(all_settings) * runs );

mean_o = mean([mean_o_high,mean_o_rest]);
std_o = mean([std_o_high,std_o_rest]);

mean_p_high = 0;
std_p_high = 0;

mean_p_rest = 0;
std_p_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_p_high = mean_p_high + mean(P(iset).run(irun).rho_high(:));
        std_p_high = std_p_high + std(P(iset).run(irun).rho_high(:));
    
        mean_p_rest = mean_p_rest + mean(P(iset).run(irun).rho_rest(:));
        std_p_rest = std_p_rest + std(P(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_p_high = mean_p_high / ( length(all_settings) * runs );
std_p_high = std_p_high / ( length(all_settings) * runs );

mean_p_rest = mean_p_rest / ( length(all_settings) * runs );
std_p_rest = std_p_rest / ( length(all_settings) * runs );

mean_p = mean([mean_p_high,mean_p_rest]);
std_p = mean([std_p_high,std_p_rest]);

mean_f_high = 0;
std_f_high = 0;

mean_f_rest = 0;
std_f_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_f_high = mean_f_high + mean(F(iset).run(irun).rho_high(:));
        std_f_high = std_f_high + std(F(iset).run(irun).rho_high(:));
    
        mean_f_rest = mean_o_rest + mean(F(iset).run(irun).rho_rest(:));
        std_f_rest = std_f_rest + std(F(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_f_high = mean_f_high / ( length(all_settings) * runs );
std_f_high = std_f_high / ( length(all_settings) * runs );

mean_f_rest = mean_f_rest / ( length(all_settings) * runs );
std_f_rest = std_f_rest / ( length(all_settings) * runs );

mean_f = mean([mean_f_high,mean_f_rest]);
std_f = mean([std_f_high,std_f_rest]);

fs = 12;

for iset=1:length(all_settings)
    
    for irun=1:runs
        
       %CLIM = [0 (mean_o + std_o)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(O(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Occipital','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(O(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Occipital','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       %CLIM = [0 (mean_p + std_p)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(P(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Parietal','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(P(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Parietal','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       %CLIM = [0 (mean_f + std_f)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(F(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Frontal','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(F(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Frontal','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
        
    end
    
    
    
end


end

function plot_FC_crosslobes(all_settings)

for iset=1:length(all_settings)
    
    settings = all_settings(iset).settings;
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

    F_n = FC_n_FO.area1_n;
    O_n = FC_n_FO.area2_n;
    P_n = FC_n_PO.area1_n;

    for irun=1:2
        
        FO(iset).run(irun).rho_high = FC_SQ_FO.run(irun).rho_high_roi(1:F_n,F_n+1:F_n+O_n);
        FO(iset).run(irun).rho_rest = FC_SQ_FO.run(irun).rho_rest_roi(1:F_n,F_n+1:F_n+O_n);
    
        FP(iset).run(irun).rho_high = FC_SQ_FP.run(irun).rho_high_roi(1:F_n,F_n+1:F_n+P_n);
        FP(iset).run(irun).rho_rest = FC_SQ_FP.run(irun).rho_rest_roi(1:F_n,F_n+1:F_n+P_n);
    
        PO(iset).run(irun).rho_high = FC_SQ_PO.run(irun).rho_high_roi(1:P_n,P_n+1:P_n+O_n);
        PO(iset).run(irun).rho_rest = FC_SQ_PO.run(irun).rho_rest_roi(1:P_n,P_n+1:P_n+O_n);
    
    end    
  
end

runs = 2;

mean_fo_high = 0;
std_fo_high = 0;

mean_fo_rest = 0;
std_fo_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_fo_high = mean_fo_high + mean(FO(iset).run(irun).rho_high(:));
        std_fo_high = std_fo_high + std(FO(iset).run(irun).rho_high(:));
    
        mean_fo_rest = mean_fo_rest + mean(FO(iset).run(irun).rho_rest(:));
        std_fo_rest = std_fo_rest + std(FO(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_fo_high = mean_fo_high / ( length(all_settings) * runs );
std_fo_high = std_fo_high / ( length(all_settings) * runs );

mean_fo_rest = mean_fo_rest / ( length(all_settings) * runs );
std_fo_rest = std_fo_rest / ( length(all_settings) * runs );

mean_fo = mean([mean_fo_high,mean_fo_rest]);
std_fo = mean([std_fo_high,std_fo_rest]);

mean_fp_high = 0;
std_fp_high = 0;

mean_fp_rest = 0;
std_fp_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_fp_high = mean_fp_high + mean(FP(iset).run(irun).rho_high(:));
        std_fp_high = std_fp_high + std(FP(iset).run(irun).rho_high(:));
    
        mean_fp_rest = mean_fp_rest + mean(FP(iset).run(irun).rho_rest(:));
        std_fp_rest = std_fp_rest + std(FP(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_fp_high = mean_fp_high / ( length(all_settings) * runs );
std_fp_high = std_fp_high / ( length(all_settings) * runs );

mean_fp_rest = mean_fp_rest / ( length(all_settings) * runs );
std_fp_rest = std_fp_rest / ( length(all_settings) * runs );

mean_fp = mean([mean_fp_high,mean_fp_rest]);
std_fp = mean([std_fp_high,std_fp_rest]);

mean_po_high = 0;
std_po_high = 0;

mean_po_rest = 0;
std_po_rest = 0;

for iset=1:length(all_settings)
    
    for irun=1:2
        
        mean_po_high = mean_po_high + mean(PO(iset).run(irun).rho_high(:));
        std_po_high = std_po_high + std(PO(iset).run(irun).rho_high(:));
    
        mean_po_rest = mean_fo_rest + mean(PO(iset).run(irun).rho_rest(:));
        std_po_rest = std_po_rest + std(PO(iset).run(irun).rho_rest(:));
    
    end
    
end

mean_po_high = mean_po_high / ( length(all_settings) * runs );
std_po_high = std_po_high / ( length(all_settings) * runs );

mean_po_rest = mean_po_rest / ( length(all_settings) * runs );
std_po_rest = std_po_rest / ( length(all_settings) * runs );

mean_po = mean([mean_po_high,mean_po_rest]);
std_po = mean([std_po_high,std_po_rest]);

fs = 12;

for iset=1:length(all_settings)
    
    for irun=1:runs
        
       %CLIM = [0 (mean_fo + std_fo)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(FO(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Frontal-Occipital','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(FO(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Frontal-Occipital','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       %CLIM = [0 (mean_fp + std_fp)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(FP(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Fronta-Parietal','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(FP(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Frontal-Parietal','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Frontal-Parietal-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       %CLIM = [0 (mean_po + std_po)];
       CLIM = [0 300];
        
       f = figure;
       
       imagesc(PO(iset).run(irun).rho_high,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Parietal-Occipital','-','High-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','HighAttention-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
       f = figure;
       
       imagesc(PO(iset).run(irun).rho_rest,CLIM);
       set(gca,'FontSize',fs);
       
       title(strcat('Parietal-Occipital','-','Rest-Run','-',int2str(irun)));
       
       print(f,'-djpeg',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.jpg'));
       print(f,'-depsc',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.eps'));
       print(f,'-dpdf',strcat(all_settings(iset).settings.folders.experiment,'-',all_settings(iset).settings.folders.subject,'-','FC-IC-SQ-AAL','-','RestingState-Parietal-Occipital-Run-',int2str(irun),'-Max-',int2str(CLIM(2)),'.pdf'));
        
        
    end
    
    
    
end

end