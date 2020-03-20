function ROI_labels = get_ROI_labels(rois)

nROI = size(rois,1);

ROI_labels = char.empty;

for iROI=1:nROI
    
    label = rois(iROI,:);
    
    label = strtrim(label);
    
    idx_start = strfind(label,'\');

    label = label(idx_start(end)+1+4:end-8);
    
    ROI_labels{iROI} = label;
    
    clear label
    
end

end
