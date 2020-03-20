
function real_FC_voxel_AAL_ROI_densities_overlaps_with_ttest


filename{1} = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both';
filename{2} = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both';
filename{3} = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both';
filename{4} = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both';

density_label{1} = 'Passive-Neg';
density_label{2} = 'Passive-Pos';
density_label{3} = 'Track-Neg';
density_label{4} = 'Track-Pos';

for iFile=1:4
    computeOverlapVolume(filename{iFile},density_label{iFile});
end

end

function computeOverlapVolume(filename,density_label)

%%% LOAD T-TEST

pcriterion = 1.6;

load_img = nifti('ttest-activation.nii');
load_img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FEAT\FEAT-ATTENTION-RESTING\',load_img.dat.fname);

TTEST_img = load_img.dat(:,:,:);

idx_attention = find(TTEST_img > pcriterion);
idx_resting = find(TTEST_img < -pcriterion);

%%% LOAD lFCD
        
load_img = nifti(strcat(filename,'.nii'));
folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\rendering\both\clusters-dlh=1.5-both';
load_img.dat.fname = strcat(folder,'/',load_img.dat.fname);

my_img = load_img.dat(:,:,:);

idx_negative = find(my_img<0);
idx_positive = find(my_img>0);

new_img = uint8(zeros(size(my_img)));

idx_attention_decrease = idx_attention(ismember(idx_attention,idx_negative));
idx_attention_increase = idx_attention(ismember(idx_attention,idx_positive));

idx_resting_decrease = idx_resting(ismember(idx_resting,idx_negative));
idx_resting_increase = idx_resting(ismember(idx_resting,idx_positive));

code = 1;
for iidx=1:length(idx_attention_decrease)

    idx = idx_attention_decrease(iidx);
    
    [x,y,z] = ind2sub(size(my_img),idx);
    
    new_img(x,y,z) = code;

end

code = 2;
for iidx=1:length(idx_attention_increase)

    idx = idx_attention_increase(iidx);
    
    [x,y,z] = ind2sub(size(my_img),idx);
    
    new_img(x,y,z) = code;

end

code = 3;
for iidx=1:length(idx_resting_decrease)

    idx = idx_resting_decrease(iidx);
    
    [x,y,z] = ind2sub(size(my_img),idx);
    
    new_img(x,y,z) = code;

end

code = 4;
for iidx=1:length(idx_resting_increase)

    idx = idx_resting_increase(iidx);
    
    [x,y,z] = ind2sub(size(my_img),idx);
    
    new_img(x,y,z) = code;

end

nifti_file = load_img;
offset = load_img.dat.offset;
scl_slope = load_img.dat.scl_slope;
scl_inter = load_img.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_img.dat.dim;

descrip = 'Overlaps';

fname = strcat('LHR','-','All-Subjects','-','Overlaps-TTest-lFCD','-',density_label,'.nii');
input_data = new_img; 
real_save_image;


end

