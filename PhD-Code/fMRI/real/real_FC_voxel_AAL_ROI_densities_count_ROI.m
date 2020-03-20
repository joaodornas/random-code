

filename{1} = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both';
filename{2} = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both';
filename{3} = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both';
filename{4} = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both';

density_label{1} = 'Passive-Neg';
density_label{2} = 'Passive-Pos';
density_label{3} = 'Track-Neg';
density_label{4} = 'Track-Pos';

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

column_decrease = -1;
column_increase = 0;
for ifile=1:4
    
    load_img = nifti(strcat(filename{ifile},'.nii'));
    folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\rendering\both\clusters-dlh=1.5-both';
    load_img.dat.fname = strcat(folder,'/',load_img.dat.fname);

    my_img = load_img.dat(:,:,:);

    idx_decrease = find(my_img<0);
    idx_increase = find(my_img>0);
    
    column_decrease = column_decrease + 2;
    column_increase = column_increase + 2;
    
    countROI{1,column_decrease+1} = strcat(density_label{ifile},'-','decrease');
    countROI{1,column_increase+1} = strcat(density_label{ifile},'-','increase');

    for iROI=1:nNodes
        
        countROI{iROI+1,1} = AAL_ROI(iROI).Nom_L;
        
        idx_ROI = AAL_ROI(iROI).ID;
        idx_voxels = find(AAL_img==idx_ROI);

        countROI{iROI+1,column_decrease+1} = length(find(ismember(idx_decrease,idx_voxels)));
        countROI{iROI+1,column_increase+1} = length(find(ismember(idx_increase,idx_voxels)));
  
    end
    
end
