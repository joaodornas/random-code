function plot_my_slices(img,MNI_img)

figure;

%%% TRANSVERSAL - z

subplot(3,3,1);

slice = squeeze(img(1:end,1:end,55));
MNI_slice = squeeze(MNI_img(1:end,1:end,65));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,2);

slice = squeeze(img(1:end,1:end,45));
MNI_slice = squeeze(MNI_img(1:end,1:end,45));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,3);

slice = squeeze(img(1:end,1:end,55));
MNI_slice = squeeze(MNI_img(1:end,1:end,25));

h = imagesc(slice);

colormap jet;
colorbar;

%%% CORONAL - x

subplot(3,3,4);

slice = squeeze(img(25,1:end,1:end));
MNI_slice = squeeze(MNI_img(25,1:end,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,5);

slice = squeeze(img(45,1:end,1:end));
MNI_slice = squeeze(MNI_img(45,1:end,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,6);

slice = squeeze(img(65,1:end,1:end));
MNI_slice = squeeze(MNI_img(65,1:end,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

%%% SAGITTAL - y

subplot(3,3,7);

slice = squeeze(img(1:end,34,1:end));
MNI_slice = squeeze(MNI_img(1:end,34,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,8);

slice = squeeze(img(1:end,54,1:end));
MNI_slice = squeeze(MNI_img(1:end,54,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

subplot(3,3,9);

slice = squeeze(img(1:end,74,1:end));
MNI_slice = squeeze(MNI_img(1:end,74,1:end));

h = imagesc(slice);

colormap jet;
colorbar;

end

