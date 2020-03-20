function lowhigh_FC_ROI

settings_jan_0502;

[rho, pval, roislabels] = doTheMath(settings);
plotFC(rho,pval,roislabels,settings);


return

function [rho, pval, roislabels] = doTheMath(settings)

get_at_this_preprocessed_step = settings.folders.normalization.name;
prefix_for_the_preprocessed_step = settings.folders.normalization.prefix;

settings = lowhigh_concatenate_folders_strings(settings,get_at_this_preprocessed_step);

roi_files_cortex_RL = get_ROI_aal_files_cortex_RL;
rois = maroi('load_cell', roi_files_cortex_RL);

roislabels = get_ROI_labels(roi_files_cortex_RL);

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

function plotFC(rho,pval,roislabels,settings)


h1 = figure;

subplot(3,2,1);
imagesc(rho.restingstate.run1);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('Resting State - Run 1');
colorbar;

subplot(3,2,2);
imagesc(rho.restingstate.run2);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('Resting State - Run 2');
colorbar;

subplot(3,2,3);
imagesc(rho.mot4.run1);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('High Attention - Run 1');
colorbar;

subplot(3,2,4);
imagesc(rho.mot4.run2);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('High Attention - Run 2');
colorbar;

subplot(3,2,5);
imagesc(rho.mot2.run1);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('Low Attention - Run 1');
colorbar;

subplot(3,2,6);
imagesc(rho.mot2.run2);
% set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
% xticklabel_rotate(1:length(roislabels),90,roislabels);
title('Low Attention - Run 2');
colorbar;

print(h1, '-depsc', strcat(settings.folders.subject,'-','FC'));
print(h1, '-djpeg', strcat(settings.folders.subject,'-','FC'));


f = figure;

imagesc(rho.restingstate.run1);
title('Resting State - Run 1');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(f, '-depsc', strcat(settings.folders.subject,'-RestingState-','FC'));
print(f, '-djpeg', strcat(settings.folders.subject,'-RestingState-','FC'));

g = figure;

imagesc(rho.restingstate.run2);
title('Resting State - Run 2');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(g, '-depsc', strcat(settings.folders.subject,'-RestingState-','FC'));
print(g, '-djpeg', strcat(settings.folders.subject,'-RestingState-','FC'));

h = figure;

imagesc(rho.mot4.run1);
title('High Attention - Run 1');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(h, '-depsc', strcat(settings.folders.subject,'-HighAttention-','FC'));
print(h, '-djpeg', strcat(settings.folders.subject,'-HighAttention-','FC'));

i = figure;

imagesc(rho.mot4.run2);
title('High Attention - Run 2');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(i, '-depsc', strcat(settings.folders.subject,'-HighAttention-','FC'));
print(i, '-djpeg', strcat(settings.folders.subject,'-HighAttention-','FC'));

j = figure;

imagesc(rho.mot2.run1);
title('Low Attention - Run 1');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(j, '-depsc', strcat(settings.folders.subject,'-LowAttention-','FC'));
print(j, '-djpeg', strcat(settings.folders.subject,'-LowAttention-','FC'));

k = figure;

imagesc(rho.mot2.run2);
title('Low Attention - Run 2');
colorbar;

set(gca,'FontSize',3,'XTickLabel',roislabels,'YTick',1:length(roislabels),'YTickLabel',roislabels');
xticklabel_rotate(1:length(roislabels),90,roislabels);

print(k, '-depsc', strcat(settings.folders.subject,'-LowAttention-','FC'));
print(k, '-djpeg', strcat(settings.folders.subject,'-LowAttention-','FC'));

return

