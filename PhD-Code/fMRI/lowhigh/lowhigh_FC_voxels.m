function lowhigh_FC_voxels

settings_jan_0502;

[rho, pval, rois] = doTheMath(settings);
plotFC(rho,pval, rois,settings);


return

function [rho, pval, rois] = doTheMath(settings)

get_at_this_preprocessed_step = settings.folders.normalization.name;
prefix_for_the_preprocessed_step = settings.folders.normalization.prefix;

settings = lowhigh_concatenate_folders_strings(settings,get_at_this_preprocessed_step);

images_files.mot4.run1 = get_images_files(settings.mot4.run(1).fullpath,prefix_for_the_preprocessed_step,settings.functional.mot4.run(1).prefix,settings.functional.mot4.firstTR,settings.functional.mot4.lastTR);
mY.mot4.run1 = get_marsy(rois{:},images_files.mot4.run1,'mean');
y.mot4.run1 = summary_data(mY.mot4.run1);
[rho.mot4.run1,pval.mot4.run1] = corr(y.mot4.run1,'type','pearson');

images_files.mot4.run2 = get_images_files(settings.mot4.run(2).fullpath,prefix_for_the_preprocessed_step,settings.functional.mot4.run(2).prefix,settings.functional.mot4.firstTR,settings.functional.mot4.lastTR);
mY.mot4.run2 = get_marsy(rois{:},images_files.mot4.run2,'mean');
y.mot4.run2 = summary_data(mY.mot4.run2);
[rho.mot4.run2,pval.mot4.run2] = corr(y.mot4.run2,'type','pearson');

images_files.mot2.run1 = get_images_files(settings.mot2.run(1).fullpath,prefix_for_the_preprocessed_step,settings.functional.mot2.run(1).prefix,settings.functional.mot2.firstTR,settings.functional.mot2.lastTR);
mY.mot2.run1 = get_marsy(rois{:},images_files.mot2.run1,'mean');
y.mot2.run1 = summary_data(mY.mot2.run1);
[rho.mot2.run1,pval.mot2.run1] = corr(y.mot2.run1,'type','pearson');

images_files.mot2.run2 = get_images_files(settings.mot2.run(2).fullpath,prefix_for_the_preprocessed_step,settings.functional.mot2.run(2).prefix,settings.functional.mot2.firstTR,settings.functional.mot2.lastTR);
mY.mot2.run2 = get_marsy(rois{:},images_files.mot2.run2,'mean');
y.mot2.run2 = summary_data(mY.mot2.run2);
[rho.mot2.run2,pval.mot2.run2] = corr(y.mot2.run2,'type','pearson');

images_files.restingstate.run1 = get_images_files(settings.restingstate.run(1).fullpath,prefix_for_the_preprocessed_step,settings.functional.restingstate.run(1).prefix,settings.functional.restingstate.firstTR,settings.functional.restingstate.lastTR);
mY.restingstate.run1 = get_marsy(rois{:},images_files.restingstate.run1,'mean');
y.restingstate.run1 = summary_data(mY.restingstate.run1);
[rho.restingstate.run1,pval.restingstate.run1] = corr(y.restingstate.run1,'type','pearson');

images_files.restingstate.run2 = get_images_files(settings.restingstate.run(2).fullpath,prefix_for_the_preprocessed_step,settings.functional.restingstate.run(2).prefix,settings.functional.restingstate.firstTR,settings.functional.restingstate.lastTR);
mY.restingstate.run2 = get_marsy(rois{:},images_files.restingstate.run2,'mean');
y.restingstate.run2 = summary_data(mY.restingstate.run2);
[rho.restingstate.run2,pval.restingstate.run2] = corr(y.restingstate.run2,'type','pearson');

return

function plotFC(rho,pval,rois,settings)

% % ROI_labels = get_ROI_labels(rois);

h = figure;

subplot(3,2,1);
imagesc(rho.restingstate.run1);
title('Resting State - Run 1');
colorbar;

subplot(3,2,2);
imagesc(rho.restingstate.run2);
title('Resting State - Run 2');
colorbar;

subplot(3,2,3);
imagesc(rho.mot4.run1);
title('High Attention - Run 1');
colorbar;

subplot(3,2,4);
imagesc(rho.mot4.run2);
title('High Attention - Run 2');
colorbar;

subplot(3,2,5);
imagesc(rho.mot2.run1);
title('Low Attention - Run 1');
colorbar;

subplot(3,2,6);
imagesc(rho.mot2.run2);
title('Low Attention - Run 2');
colorbar;

print(h, '-depsc', strcat(settings.folders.subject,'-','FC'));
print(h, '-djpeg', strcat(settings.folders.subject,'-','FC'));

return

