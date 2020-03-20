function teste(MOT4_Threshold,MOT2_Threshold,RestingState_Threshold,MOT4RestingStateDiff,MOT2RestingStateDiff,MOT4MOT2Diff,pvalue,kind,common_mask,common_idx,settings,scale)

if strcmpi('threshold',kind)
    
    %mask_vec = reshape(common_mask,[1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3)]);
    %all_idx = find(common_mask<2);
    %all_idx(common_idx) = [];
    idx_brain = common_idx;
    
else
    
    %mask_vec = reshape(common_mask,[1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3)]);
    idx_brain = find(common_mask ~= 0);

end

%MOT4_threshold_vec = reshape(MOT4_Threshold,[1,size(MOT4_Threshold,1)*size(MOT4_Threshold,2)*size(MOT4_Threshold,3)]);
% MOT4_threshold_vec(idx_space) = [];
MOT4_threshold_vec = double(MOT4_Threshold(idx_brain))*10^(-scale);
%
%MOT2_threshold_vec = reshape(MOT2_Threshold,[1,size(MOT2_Threshold,1)*size(MOT2_Threshold,2)*size(MOT2_Threshold,3)]);
% MOT2_threshold_vec(idx_space) = [];
MOT2_threshold_vec = double(MOT2_Threshold(idx_brain))*10^(-scale);
% 
%RestingState_threshold_vec = reshape(RestingState_Threshold,[1,size(RestingState_Threshold,1)*size(RestingState_Threshold,2)*size(RestingState_Threshold,3)]);
% RestingState_threshold_vec(idx_space) = [];
RestingState_threshold_vec = double(RestingState_Threshold(idx_brain))*10^(-scale);
% 
%MOT4RestingStateDiff_vec = reshape(MOT4RestingStateDiff,[1,size(MOT4RestingStateDiff,1)*size(MOT4RestingStateDiff,2)*size(MOT4RestingStateDiff,3)]);
% MOT4RestingStateDiff_vec(idx_space) = [];
MOT4RestingStateDiff_vec = double(MOT4RestingStateDiff(idx_brain))*10^(-scale);
% 
%MOT2RestingStateDiff_vec = reshape(MOT2RestingStateDiff,[1,size(MOT2RestingStateDiff,1)*size(MOT2RestingStateDiff,2)*size(MOT2RestingStateDiff,3)]);
% MOT2RestingStateDiff_vec(idx_space) = [];
MOT2RestingStateDiff_vec = double(MOT2RestingStateDiff(idx_brain))*10^(-scale);
% 
%MOT4MOT2Diff_vec = reshape(MOT4MOT2Diff,[1,size(MOT4MOT2Diff,1)*size(MOT4MOT2Diff,2)*size(MOT4MOT2Diff,3)]);
% MOT4MOT2Diff_vec(idx_space) = [];
MOT4MOT2Diff_vec = double(MOT4MOT2Diff(idx_brain))*10^(-scale);

MOT4_RestingState_covariance = getCovariance(MOT4_threshold_vec,RestingState_threshold_vec);
MOT2_RestingState_covariance = getCovariance(MOT2_threshold_vec,RestingState_threshold_vec);
MOT4_MOT2_covariance = getCovariance(MOT4_threshold_vec,MOT2_threshold_vec);

max_threshold = ceil(max([max(MOT4_threshold_vec),max(MOT2_threshold_vec),max(RestingState_threshold_vec)]));
min_threshold = floor(min([min(MOT4_threshold_vec),min(MOT2_threshold_vec),min(RestingState_threshold_vec)]));

max_diff = ceil(max([max(MOT4RestingStateDiff_vec),max(MOT2RestingStateDiff_vec),max(MOT4MOT2Diff_vec)]));
min_diff = floor(min([min(MOT4RestingStateDiff_vec),min(MOT2RestingStateDiff_vec),min(MOT4MOT2Diff_vec)]));

nPoints = 1000;

bins_centers_threshold = linspace(min_threshold,max_threshold,nPoints);
bins_centers_diff = linspace(min_diff,max_diff,nPoints);

f = figure;
hold on;
fs = 6;

%%% PLOT Z-SCORE THRESHOLDED

[nquality,xquality] = hist(MOT4_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_MOT4 = fquality;
xquality_MOT4 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4_threshold_vec = double(mean(MOT4_Threshold(idx_brain),'native'))*10^(-scale);
X_MOT4 = mean_x_MOT4_threshold_vec;
Y_MOT4 = max(dquality_MOT4) / 2;

[nquality,xquality] = hist(MOT2_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_MOT2 = fquality;
xquality_MOT2 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT2_threshold_vec = double(mean(MOT2_Threshold(idx_brain),'native'))*10^(-scale);
X_MOT2 = mean_x_MOT2_threshold_vec;
Y_MOT2 = max(dquality_MOT2) / 2;

[nquality,xquality] = hist(RestingState_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_RestingState = fquality;
xquality_RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_RestingState_threshold_vec = double(mean(RestingState_Threshold(idx_brain),'native'))*10^(-scale);
X_RestingState = mean_x_RestingState_threshold_vec;
Y_RestingState = max(dquality_RestingState) / 2;

min_dquality = min([dquality_MOT4, dquality_MOT2, dquality_RestingState]);
max_dquality = max([dquality_MOT4, dquality_MOT2, dquality_RestingState]);

subplot(2,3,1);
bar( xquality_MOT4, dquality_MOT4  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(X_MOT4,Y_MOT4,strcat('mean:',num2str(mean_x_MOT4_threshold_vec)));
title('MOT 4');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

subplot(2,3,2);
bar( xquality_MOT2, dquality_MOT2  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(X_MOT2,Y_MOT2,strcat('mean:',num2str(mean_x_MOT2_threshold_vec)));
title('MOT 2');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

subplot(2,3,3);
bar( xquality_RestingState, dquality_RestingState  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(X_RestingState,Y_RestingState,strcat('mean:',num2str(mean_x_RestingState_threshold_vec)));
title('Resting State');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

%%% PLOT DIFFERENCE

[nquality,xquality] = hist(MOT4RestingStateDiff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT4RestingState = fquality;
xquality_MOT4RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4RestingState = double(mean(MOT4RestingStateDiff(idx_brain),'native'))*10^(-scale);
X_MOT4RestingState = mean_x_MOT4RestingState;
Y_MOT4RestingState = max(dquality_MOT4RestingState) / 2;

[nquality,xquality] = hist(MOT2RestingStateDiff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT2RestingState = fquality;
xquality_MOT2RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT2RestingState = double(mean(MOT2RestingStateDiff(idx_brain),'native'))*10^(-scale);
X_MOT2RestingState = mean_x_MOT2RestingState;
Y_MOT2RestingState = max(dquality_MOT2RestingState) / 2;

[nquality,xquality] = hist(MOT4MOT2Diff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT4MOT2 = fquality;
xquality_MOT4MOT2 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4MOT2 = double(mean(MOT4MOT2Diff(idx_brain),'native'))*10^(-scale);
X_MOT4MOT2 = mean_x_MOT4MOT2;
Y_MOT4MOT2 = max(dquality_MOT4MOT2) / 2;

min_dquality = min([dquality_MOT4RestingState, dquality_MOT2RestingState, dquality_MOT4MOT2]);
max_dquality = max([dquality_MOT4RestingState, dquality_MOT2RestingState, dquality_MOT4MOT2]);

subplot(2,3,4);
bar( xquality_MOT4RestingState, dquality_MOT4RestingState  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 4 - Resting State');
text(X_MOT4RestingState,Y_MOT4RestingState,strcat('mean:',num2str(mean_x_MOT4RestingState)));
text(X_MOT4RestingState,Y_MOT4RestingState/2,strcat('cov:',num2str(MOT4_RestingState_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

subplot(2,3,5);
bar( xquality_MOT2RestingState, dquality_MOT2RestingState  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 2 - Resting State');
text(X_MOT2RestingState,Y_MOT2RestingState,strcat('mean:',num2str(mean_x_MOT2RestingState)));
text(X_MOT2RestingState,Y_MOT2RestingState/2,strcat('cov:',num2str(MOT2_RestingState_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

subplot(2,3,6);
bar( xquality_MOT4MOT2, dquality_MOT4MOT2  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 4 - MOT 2');
text(X_MOT4MOT2,Y_MOT4MOT2,strcat('mean:',num2str(mean_x_MOT4MOT2)));
text(X_MOT4MOT2,Y_MOT4MOT2/2,strcat('cov:',num2str(MOT4_MOT2_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

if strcmpi('threshold',kind)
    
    suptitle(strcat('p-value:',num2str(pvalue)));

else
    
    suptitle('No Threshold');
    
end

print(f,'-djpeg',strcat('Low-High-',settings.folders.subject,'-','ZScore-Diff-Histogram-AllRuns-p-value-',num2str(pvalue),'.jpeg'));
print(f,'-depsc',strcat('Low-High-',settings.folders.subject,'-','ZScore-Diff-Histogram-AllRuns-p-value-',num2str(pvalue),'.eps'));

function covariance = getCovariance(X,Y)

diagonal = diag(flip(cov(X,Y)));
covariance = diagonal(1);

end

end