
settings_jan_0805;

get_at_this_preprocessed_step = settings.FSL.folders.custom;

file = settings.FSL.files.functional.custom.filtered;

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

run = 2;
filename = strcat(settings.mot4.FSL.run(run).fullpath,'\',file);

datafile = nifti(strcat(filename,'.nii'));
datafile.dat.fname = strcat(filename,'.nii');
data = datafile.dat(:,:,:,:);

nTR = size(data,4);

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

iNode = 4;

idx_ROI = AAL_ROI(idx_nodes(iNode)).ID;

idx_voxels = find(AAL_img == idx_ROI);

nVoxels = length(idx_voxels);

time_series = zeros(nVoxels,nTR);
for iVoxel=1:nVoxels

    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

    time_series(iVoxel,:) = squeeze(data(idxx,idxy,idxz,:));

end

%plot(time_series);
%xlim([0 nTR]);

%threshold = [3 5 10 15 20];
threshold = [20 30 50 70 90].*3;
wavelet = 'd4';
nscale = {'conservative','liberal','extreme'};

for i=1:5
    
    for iscale=2
 
        clean(i).ts = zeros(nVoxels,nTR);
        noise(i).ts = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
        
            [clean(i).ts(iVoxel,:), noise(i).ts(iVoxel,:)] = WaveletDespike(squeeze(time_series(iVoxel,:)),'wavelet',wavelet,'threshold',threshold(i),'nscale',nscale{iscale});
    
        end
        
    end
    
end

max_value = max(mean(time_series,1));
min_value = min(mean(time_series,1));

figure;

k = 0;

for i=1:5
    
    for iscale=2
        
        k = k + 1;
        
        subplot(5,3,k);
        
        plot(mean(time_series,1),'r');
        xlim([0 nTR]);
        ylim([min_value max_value]);
        ylabel('BOLD');
        xlabel('TR');
        legend({'Original'});
        title(strcat('scale:',nscale{iscale},'-','threshold:',int2str(threshold(i))));
        
        k = k + 1;
        
        subplot(5,3,k);
        
        plot(mean(clean(i).ts,1),'b');
        xlim([0 nTR]);
        ylim([min_value max_value]);
        ylabel('BOLD');
        xlabel('TR');
        legend({'Cleaned'});
        title(strcat('scale:',nscale{iscale},'-','threshold:',int2str(threshold(i))));
        
        k = k + 1;
        
        subplot(5,3,k);
        
        plot(mean(time_series,1),'r');
        hold on;
        plot(mean(clean(i).ts,1),'b');
        xlim([0 nTR]);
        ylim([min_value max_value]);
        ylabel('BOLD');
        xlabel('TR');
        legend({'Original','Cleaned'});
        title(strcat('scale:',nscale{iscale},'-','threshold:',int2str(threshold(i))));
        
    end
    
end
