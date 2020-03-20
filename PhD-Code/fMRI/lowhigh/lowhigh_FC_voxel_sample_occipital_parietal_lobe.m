function lowhigh_FC_voxel_sample_occipital_parietal_lobe


settings_jan_0805;
%settings_elena_2905;

doTheMath(settings);

%plotResults(settings);

%plotCumulative(settings);

end

function doTheMath(settings)

file = settings.FSL.files.functional.custom.residual_voxel;

get_at_this_preprocessed_step = settings.FSL.folders.custom;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nTR = size(MOT4Run1,4);

idx_occipital_lobe = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal_lobe = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

Voxels_label{1} = 'High-Attention-Run-1';
Voxels_label{2} = 'High-Attention-Run-2';
Voxels_label{3} = 'Low-Attention-Run-1';
Voxels_label{4} = 'Low-Attention-Run-2';
Voxels_label{5} = 'RestingState-Run-1';
Voxels_label{6} = 'RestingState-Run-2';

all_idx_voxels_occipital = [];
for iidx=1:length(idx_occipital_lobe)
    
   idx = idx_occipital_lobe(iidx);
   idx_ROI = AAL_ROI(idx).ID;
   
   idx_voxels = find(AAL_img==idx_ROI);
   if size(idx_voxels,1) ~= 1, idx_voxels = idx_voxels'; end
   
   all_idx_voxels_occipital = [all_idx_voxels_occipital, idx_voxels];
   
end

all_idx_voxels_parietal = [];
for iidx=1:length(idx_parietal_lobe)
    
   idx = idx_parietal_lobe(iidx);
   idx_ROI = AAL_ROI(idx).ID;
   
   idx_voxels = find(AAL_img==idx_ROI);
   if size(idx_voxels,1) ~= 1, idx_voxels = idx_voxels'; end
   
   all_idx_voxels_parietal = [all_idx_voxels_parietal, idx_voxels];
   
end

nOccipital = length(all_idx_voxels_occipital);
nParietal = length(all_idx_voxels_parietal);

nVoxelsSamples = 1000;
nComputations = 10;

for isample=1:nComputations
    
    Voxels_samples(isample).idx_occipital = randi(nOccipital,nVoxelsSamples,1);
    Voxels_samples(isample).idx_parietal = randi(nParietal,nVoxelsSamples,1);
        
end

for isample=1:nComputations

    MOT4Run1_lobe = zeros(nVoxelsSamples*2,nTR);
    MOT4Run2_lobe = zeros(nVoxelsSamples*2,nTR);
    MOT2Run1_lobe = zeros(nVoxelsSamples*2,nTR);
    MOT2Run2_lobe = zeros(nOccipital+nParietal,nTR);
    RestingStateRun1_lobe = zeros(nVoxelsSamples*2,nTR);
    RestingStateRun2_lobe = zeros(nVoxelsSamples*2,nTR);

    for iidx=1:nVoxelsSamples
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels_occipital(Voxels_samples(isample).idx_occipital(iidx)));
   
        MOT4Run1_lobe(iidx,:) = MOT4Run1(idxx,idxy,idxz,:);
        MOT4Run2_lobe(iidx,:) = MOT4Run2(idxx,idxy,idxz,:);
        MOT2Run1_lobe(iidx,:) = MOT2Run1(idxx,idxy,idxz,:);
        MOT2Run2_lobe(iidx,:) = MOT2Run2(idxx,idxy,idxz,:);
        RestingStateRun1_lobe(iidx,:) = RestingStateRun1(idxx,idxy,idxz,:);
        RestingStateRun2_lobe(iidx,:) = RestingStateRun2(idxx,idxy,idxz,:);
    
    end

    for iidx=1:nVoxelsSamples
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels_parietal(Voxels_samples(isample).idx_parietal(iidx)));
   
        MOT4Run1_lobe(nVoxelsSamples+iidx,:) = MOT4Run1(idxx,idxy,idxz,:);
        MOT4Run2_lobe(nVoxelsSamples+iidx,:) = MOT4Run2(idxx,idxy,idxz,:);
        MOT2Run1_lobe(nVoxelsSamples+iidx,:) = MOT2Run1(idxx,idxy,idxz,:);
        MOT2Run2_lobe(nVoxelsSamples+iidx,:) = MOT2Run2(idxx,idxy,idxz,:);
        RestingStateRun1_lobe(nVoxelsSamples+iidx,:) = RestingStateRun1(idxx,idxy,idxz,:);
        RestingStateRun2_lobe(nVoxelsSamples+iidx,:) = RestingStateRun2(idxx,idxy,idxz,:);
    
    end

    Voxels_nSamples = nVoxelsSamples;

    Voxels_FC(1).samples(isample).corr = gpu_correlation(MOT4Run1_lobe,Voxels_label{1},nVoxelsSamples*2);
    Voxels_FC(2).samples(isample).corr = gpu_correlation(MOT4Run2_lobe,Voxels_label{2},nVoxelsSamples*2);
    Voxels_FC(3).samples(isample).corr = gpu_correlation(MOT2Run1_lobe,Voxels_label{3},nVoxelsSamples*2);
    Voxels_FC(4).samples(isample).corr = gpu_correlation(MOT2Run2_lobe,Voxels_label{4},nVoxelsSamples*2);
    Voxels_FC(5).samples(isample).corr = gpu_correlation(RestingStateRun1_lobe,Voxels_label{5},nVoxelsSamples*2);
    Voxels_FC(6).samples(isample).corr = gpu_correlation(RestingStateRun2_lobe,Voxels_label{6},nVoxelsSamples*2);
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Voxel-Occipital-Parietal-lobe-All-Runs','.mat'),'Voxels_label','Voxels_nSamples','Voxels_FC','Voxels_samples');

end
    
function corr_mat = gpu_correlation(RUN,label,nVoxels)
        
disp(label);

corr_mat = zeros(nVoxels,nVoxels);
%pval_mat = zeros(nVoxels,nVoxels);

parfor_progress(nVoxels);

parfor iVoxel=1:nVoxels
    
   %tic
    
   X = RUN(iVoxel,:);
   gpu_X = gpuArray(X);
   if size(gpu_X,1) == 1; gpu_X = gpu_X'; end
   
   for iiVoxel=1:nVoxels
      
       Y = RUN(iiVoxel,:);
       gpu_Y = gpuArray(Y);
       if size(gpu_Y,1) == 1; gpu_Y = gpu_Y'; end
       
       %[gpu_rho,gpu_pval] = corrcoef(gpu_X,gpu_Y);
       gpu_rho = corrcoef(gpu_X,gpu_Y);
       rho = gather(gpu_rho(1,2));
       %pval = gather(gpu_pval(1,2));
       
       corr_mat(iVoxel,iiVoxel) = rho;
       %pval_mat(iVoxel,iiVoxel) = pval;
       
   end
   
   %toc
   
   parfor_progress;
    
end

end

function plotResults(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Voxel-Occipital-Parietal-lobe-All-Runs','.mat'));

for ifolder=1:length(ICA_folder)
    
   max_rho(ifolder) = max(max(ICA_FC(ifolder).rho));
   min_rho(ifolder) = min(min(ICA_FC(ifolder).rho));
    
end

max_all = max(max_rho);
min_all = min(min_rho);

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
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.pdf'));
   
   
%    f = figure;
%    
%    imagesc(ICA_FC(ifolder).percent_rho,[min_percent_all max_percent_all]);
%    
%    xlabel('occipital - parietal');
%    ylabel('parietal - occipital');
%    
%    title(strcat(ICA_label{ifolder},':',int2str(round(ICA_n(ifolder).variance_explained*100)),'% of variance'));
%    
%    colorbar;
%    
%    hold on;
%    
%    p1(1) = ICA_n(ifolder).nICpercentOccipital;
%    p2(1) = 1;
%    p1(2) = ICA_n(ifolder).nICpercentOccipital;
%    p2(2) = ICA_n(ifolder).nICpercentOccipital + ICA_n(ifolder).nICpercentParietal;
%    
%    plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
%    
%    p1(1) = 1;
%    p2(1) = ICA_n(ifolder).nICpercentOccipital;
%    p1(2) = ICA_n(ifolder).nICpercentOccipital + ICA_n(ifolder).nICpercentParietal;
%    p2(2) = ICA_n(ifolder).nICpercentOccipital;
%    
%    plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
%    
%    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.jpeg'));
%    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.eps'));
%    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.pdf'));
   
end

end

function plotCumulative(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Voxel-Occipital-Parietal-lobe-All-Runs','.mat'));

for ifolder=1:6
    
    for isample=1:Voxels_nSamples
        
        nOccipital = length(Voxels_samples(isample).idx_occipital);
        nParietal = length(Voxels_samples(isample).idx_parietal);
        
        NoA = nOccipital;
        NpA = nParietal;
   
        OOA = Voxels_FC(ifolder).samples(isample).corr(1:NoA,1:NoA);

        PPA = Voxels_FC(ifolder).samples(isample).corr(NoA+1:NoA+NpA,NoA+1:NoA+NpA);

        OPA = Voxels_FC(ifolder).samples(isample).corr(1:NoA, NoA+1:NoA+NpA); 
    
        NA = min( [NoA NpA] );  

        SOOA = zeros(isample,NA);
        SPPA = zeros(isample,NA);
        SOPA = zeros(isample,NA);

        for nsc = 2:NA    % number of components to be summed

            for ic = 2 : nsc   % loop over components in i-direction
                for jc = 1 : ic-1   % loop over componnets in j-direction

                    SOOA(isample,nsc) = SOOA(isample,nsc) + abs( OOA(ic,jc) );
                    SPPA(isample,nsc) = SPPA(isample,nsc) + abs( PPA(ic,jc) );
                    SOPA(isample,nsc) = SOPA(isample,nsc) + abs( OPA(ic,jc) );

                end
            end

            SOOA(isample,nsc) = 2 * SOOA(isample,nsc) / (nsc*nsc);   % normalization
            SPPA(isample,nsc) = 2 * SPPA(isample,nsc) / (nsc*nsc);   % normalization
            SOPA(isample,nsc) = 2 * SOPA(isample,nsc) / (nsc*nsc);   % normalization
        
        end
        
    end
    
    SOOA_mean = mean(SOOA,1);
    SPPA_mean = mean(SPPA,1);
    SOPA_mean = mean(SOPA,1);
    
    SOOA_std = std(SOOA,0,1);
    SPPA_std = std(SPPA,0,1);
    SOPA_std = std(SOPA,0,1);
    
    f = figure;
    
    fs = 14;
    
    hold on;
    shadedErrorBar(1:NA, SOOA_mean, SOOA_std, 'k');
    shadedErrorBar(1:NA, SPPA_mean, SPPA_std, 'b');
    shadedErrorBar(1:NA, SOPA_mean, SOPA_std, 'r');
    hold off;
    
    h = legend('O-O', 'P-P', 'O-P');
    set(h,'FontSize', fs, 'Location', 'NorthEast' );
    axis 'square';
    axis([0 NA 0 0.5]);
    title(ICA_label{ifolder}, 'FontSize', fs);
    xlabel( 'component no', 'FontSize', fs);
    ylabel( 'mean abs corr', 'FontSize', fs );

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Cumulative-ICA-Occipital-Parietal-lobe','-',ICA_label{ifolder},'.pdf'));
   
end
   
end

