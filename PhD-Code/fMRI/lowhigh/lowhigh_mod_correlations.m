
function lowhigh_mod_correlations

settings_jan_0805;
%settings_elena_2905;

doTheMath(settings);

%plotResults(settings);


end


function doTheMath(settings)

%%% LOAD AAL TIME SERIES

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
get_at_this_preprocessed_step = settings.FSL.folders.custom;
%file = settings.FSL.files.functional.custom.filtered;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%% LOAD DATA

lowhigh_load_all_data_FSL;

lowhigh_AAL_parcellation;

%%% LOAD MOMENTOS OF DIFFICULTY

log = load(strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\','log-settings','\','log-',settings.folders.subject,'\','high_realfmri','.mat'));

velocity = log.Block{1,1}.Velocities;
targets = log.Block{1,1}.TrialsTargets;
TrialsMOTidx = log.Block{1,1}.TrialsMOTidx;

%areaAnalyzeFolder = 'K:\Dropbox (Uni Magdeburg)\_STIMULI\LOW-HIGH-ATTENTION\PLAN-A\v5.0-29-05-2015\MOT\MOT-singleTrialBlock\AnalysesArea';
AnalyzeFolder = 'K:\Dropbox (Uni Magdeburg)\_STIMULI\LOW-HIGH-ATTENTION\PLAN-A\v5.0-29-05-2015\MOT\MOT-singleTrialBlock\Analyses';

analyze = load(strcat(AnalyzeFolder,'\','MOT-Analyze-V2-11-m-',int2str(velocity),'-',int2str(TrialsMOTidx),'.mat'));

frequency = 60;
TRinSeconds = 2;

colorEntropy = analyze.Ecolor;
targetsMotion = analyze.Wmotion(targets,:);

newColorEntropy = repmat(colorEntropy,1,analyze.hWindow);
newColorEntropy = newColorEntropy';
colorEntropy = reshape(newColorEntropy,[1,size(newColorEntropy,1)*size(newColorEntropy,2)]);

newTargetsMotion = zeros(4,length(targetsMotion)*analyze.hWindow);

for iTarget=1:4
    
    thisTarget = targetsMotion(iTarget,:);
    thisTarget = thisTarget';
    thisTarget = repmat(thisTarget,1,analyze.hWindow);
    thisTarget = thisTarget';
    thisTarget = reshape(thisTarget,[1,size(thisTarget,1)*size(thisTarget,2)]);
    
    newTargetsMotion(iTarget,:) = thisTarget(:);
    
end

colorEntropy = colorEntropy(1:frequency*TRinSeconds:end);
targetsMotion = newTargetsMotion(1:end,1:frequency*TRinSeconds:end);

%%% GET CORRELATIONS

for iCondition=1:3

    for iNode=1:nNodes 
        
        if iCondition == 1, X = mean_area_MOT4_voxel(iNode,:); end
        if iCondition == 2, X = mean_area_MOT2_voxel(iNode,:); end
        if iCondition == 3, X = mean_area_RestingState_voxel(iNode,:); end
        
        z_X = nanzscore(X);
        z_X = z_X(2:end);

        Y_color = colorEntropy(2:nTR);
        z_Y_color = nanzscore(Y_color);

        [Color_rho(iCondition,iNode),Color_pval(iCondition,iNode)] = corr(z_X',z_Y_color');

        Y_target_one = targetsMotion(1,2:nTR);    
        z_Y_target_one = nanzscore(Y_target_one);
        [Motion_one_rho(iCondition,iNode),Motion_one_pval(iCondition,iNode)] = corr(z_X',z_Y_target_one');

        Y_target_two = targetsMotion(2,2:nTR);  
        z_Y_target_two = nanzscore(Y_target_two);
        [Motion_two_rho(iCondition,iNode),Motion_two_pval(iCondition,iNode)] = corr(z_X',z_Y_target_two');

        Y_target_three = targetsMotion(3,2:nTR);  
        z_Y_target_three = nanzscore(Y_target_three);
        [Motion_three_rho(iCondition,iNode),Motion_three_pval(iCondition,iNode)] = corr(z_X',z_Y_target_three');

        Y_target_four = targetsMotion(4,2:nTR);  
        z_Y_target_four = nanzscore(Y_target_four);
        [Motion_four_rho(iCondition,iNode),Motion_four_pval(iCondition,iNode)] = corr(z_X',z_Y_target_four');

    end

end

end


function plotResults(settings)




end

function zscored = nanzscore(timeseries)

    zscored = ( timeseries - mean(timeseries(~isnan(timeseries))) ) ./ std(timeseries(~isnan(timeseries)));

end