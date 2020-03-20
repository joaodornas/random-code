function [ ROIs, labels ] = get_ROI_aal( img, summarization )

AAL_img = nifti('ROI_MNI_V4.nii');
AAL_dat = AAL_img.dat;
AAL_dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',AAL_dat.fname);
AAL_matrix = AAL_dat(:,:,:);

AAL_labels = load('ROI_MNI_V4_List.mat');

nROI = length(AAL_labels.ROI);
nTR = size(img,4);

labels = char.empty;
ROIs = zeros(nTR,nROI);

for iROI=1:nROI
   
    RGB = AAL_labels.ROI(iROI).ID;
    labels{iROI} = AAL_labels.ROI(iROI).Nom_L;
    
    idx = find(AAL_matrix == RGB);
    
    for iTR=1:nTR
        
        whole_brain = squeeze(img(:,:,:,iTR));
        
        voxels = whole_brain(idx);
        
        switch summarization
            
            case 'mean'
                
                ROI_voxels = mean(voxels);
                
        end
        
        ROIs(iTR,iROI) = ROI_voxels;
        
    end
    
end

end

