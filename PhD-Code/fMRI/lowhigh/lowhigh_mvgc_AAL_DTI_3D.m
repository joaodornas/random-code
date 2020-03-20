
function lowhigh_mvgc_AAL_DTI_3D

settings_jan_0805;
%settings_elena_2905;

filename = 'residual';

doTheMath(settings,filename);

%plotResults(settings,filename);


end

function doTheMath(settings,filename)

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
get_at_this_preprocessed_step = settings.FSL.folders.custom;
%file = settings.FSL.files.functional.custom.filtered;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%% LOAD DATA

lowhigh_load_all_data_FSL;

lowhigh_AAL_parcellation;

idx_frontal = 3;
idx_ROI = AAL_ROI(idx_frontal).ID;
idx_voxels = find(AAL_img == idx_ROI);

nVoxels = length(idx_voxels);

disp(strcat('nVoxels:',int2str(nVoxels)));

area_MOT4_voxel = zeros(nVoxels,nTR,2);
area_MOT2_voxel = zeros(nVoxels,nTR,2);
area_RestingState_voxel = zeros(nVoxels,nTR,2);

for iVoxel=1:nVoxels

    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
    processed_voxel  = preprocessVoxel(squeeze(MOT4Run1(idxx,idxy,idxz,:)));
    area_MOT4_voxel(iVoxel,:,1) = processed_voxel;

    processed_voxel  = preprocessVoxel(squeeze(MOT2Run1(idxx,idxy,idxz,:)));
    area_MOT2_voxel(iVoxel,:,1) = processed_voxel;

    processed_voxel = preprocessVoxel(squeeze(RestingStateRun1(idxx,idxy,idxz,:)));
    area_RestingState_voxel(iVoxel,:,1) = processed_voxel;

end

for iVoxel=1:nVoxels

    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

    processed_voxel  = preprocessVoxel(squeeze(MOT4Run2(idxx,idxy,idxz,:)));
    area_MOT4_voxel(iVoxel,:,2) = processed_voxel;

    processed_voxel  = preprocessVoxel(squeeze(MOT2Run2(idxx,idxy,idxz,:)));
    area_MOT2_voxel(iVoxel,:,2) = processed_voxel;

    processed_voxel = preprocessVoxel(squeeze(RestingStateRun2(idxx,idxy,idxz,:)));
    area_RestingState_voxel(iVoxel,:,2) = processed_voxel;

end

%%% GET CAUSALITY

% nTrials = [5 10 20];
% 
% for iTrial=1:3
% 
%     disp(strcat('doing trial:',int2str(nTrials(iTrial))));
%     
%     disp('Causality - High Attention');
%     reshaped_MOT4 = myReshapeRun(mean_area_MOT4_voxel,nTrials(iTrial));
%     [Granger_MOT4(iTrial).AIC,Granger_MOT4(iTrial).BIC,Granger_MOT4(iTrial).moAIC,Granger_MOT4(iTrial).moBIC,Granger_MOT4(iTrial).A,Granger_MOT4(iTrial).SIG,Granger_MOT4(iTrial).G,Granger_MOT4(iTrial).info,Granger_MOT4(iTrial).F,Granger_MOT4(iTrial).pval,Granger_MOT4(iTrial).sig] = getCausality(reshaped_MOT4,nTR,nTrials(iTrial));
% 
%     disp('Causality - Low Attention');
%     reshaped_MOT2 = myReshapeRun(mean_area_MOT2_voxel,nTrials(iTrial));
%     [Granger_MOT2(iTrial).AIC,Granger_MOT2(iTrial).BIC,Granger_MOT2(iTrial).moAIC,Granger_MOT2(iTrial).moBIC,Granger_MOT2(iTrial).A,Granger_MOT2(iTrial).SIG,Granger_MOT2(iTrial).G,Granger_MOT2(iTrial).info,Granger_MOT2(iTrial).F,Granger_MOT2(iTrial).pval,Granger_MOT2(iTrial).sig] = getCausality(reshaped_MOT2,nTR,nTrials(iTrial));
% 
%     disp('Causality - Resting State');
%     reshaped_RestingState = myReshapeRun(mean_area_RestingState_voxel,nTrials(iTrial));
%     [Granger_RestingState(iTrial).AIC,Granger_RestingState(iTrial).BIC,Granger_RestingState(iTrial).moAIC,Granger_RestingState(iTrial).moBIC,Granger_RestingState(iTrial).A,Granger_RestingState(iTrial).SIG,Granger_RestingState(iTrial).G,Granger_RestingState(iTrial).info,Granger_RestingState(iTrial).F,Granger_RestingState(iTrial).pval,Granger_RestingState(iTrial).sig] = getCausality(reshaped_RestingState,nTR,nTrials(iTrial));
% 
% end

% nTrials = 2;
% 
% disp('Causality - High Attention');
% reshaped_MOT4 = myReshapeRun(mean_area_MOT4_voxel,nTrials);
% [Granger_MOT4.AIC,Granger_MOT4.BIC,Granger_MOT4.moAIC,Granger_MOT4.moBIC,Granger_MOT4.A,Granger_MOT4.SIG,Granger_MOT4.G,Granger_MOT4.info,Granger_MOT4.F,Granger_MOT4.pval,Granger_MOT4.sig] = getCausality(reshaped_MOT4,nTR,nTrials);
% 
% disp('Causality - Low Attention');
% reshaped_MOT2 = myReshapeRun(mean_area_MOT2_voxel,nTrials);
% [Granger_MOT2.AIC,Granger_MOT2.BIC,Granger_MOT2.moAIC,Granger_MOT2.moBIC,Granger_MOT2.A,Granger_MOT2.SIG,Granger_MOT2.G,Granger_MOT2.info,Granger_MOT2.F,Granger_MOT2.pval,Granger_MOT2.sig] = getCausality(reshaped_MOT2,nTR,nTrials);
% 
% disp('Causality - Resting State');
% reshaped_RestingState = myReshapeRun(mean_area_RestingState_voxel,nTrials);
% [Granger_RestingState.AIC,Granger_RestingState.BIC,Granger_RestingState.moAIC,Granger_RestingState.moBIC,Granger_RestingState.A,Granger_RestingState.SIG,Granger_RestingState.G,Granger_RestingState.info,Granger_RestingState.F,Granger_RestingState.pval,Granger_RestingState.sig] = getCausality(reshaped_RestingState,nTR,nTrials);

% nTrials = 1;
% 
% disp('Causality - High Attention');
% reshaped_MOT4 = mean_area_MOT4_voxel;
% [Granger_MOT4.AIC,Granger_MOT4.BIC,Granger_MOT4.moAIC,Granger_MOT4.moBIC,Granger_MOT4.A,Granger_MOT4.SIG,Granger_MOT4.G,Granger_MOT4.info,Granger_MOT4.F,Granger_MOT4.pval,Granger_MOT4.sig] = getCausality(reshaped_MOT4,nTR,nTrials);
% 
% disp('Causality - Low Attention');
% reshaped_MOT2 = mean_area_MOT2_voxel;
% [Granger_MOT2.AIC,Granger_MOT2.BIC,Granger_MOT2.moAIC,Granger_MOT2.moBIC,Granger_MOT2.A,Granger_MOT2.SIG,Granger_MOT2.G,Granger_MOT2.info,Granger_MOT2.F,Granger_MOT2.pval,Granger_MOT2.sig] = getCausality(reshaped_MOT2,nTR,nTrials);
% 
% disp('Causality - Resting State');
% reshaped_RestingState = mean_area_RestingState_voxel;
% [Granger_RestingState.AIC,Granger_RestingState.BIC,Granger_RestingState.moAIC,Granger_RestingState.moBIC,Granger_RestingState.A,Granger_RestingState.SIG,Granger_RestingState.G,Granger_RestingState.info,Granger_RestingState.F,Granger_RestingState.pval,Granger_RestingState.sig] = getCausality(reshaped_RestingState,nTR,nTrials);

nTrials = 2;

for iVoxel=1:nVoxels
    
    for iiVoxel=iVoxel:nVoxels
        
        disp('Causality - High Attention');
        reshaped_MOT4(1,:,:) = area_MOT4_voxel(iVoxel,:,:);
        reshaped_MOT4(2,:,:) = area_MOT4_voxel(iiVoxel,:,:);
        
        [Granger_MOT4(iVoxel,iiVoxel).A,Granger_MOT4(iVoxel,iiVoxel).SIG,Granger_MOT4(iVoxel,iiVoxel).G,Granger_MOT4(iVoxel,iiVoxel).info,Granger_MOT4(iVoxel,iiVoxel).F,Granger_MOT4(iVoxel,iiVoxel).pval,Granger_MOT4(iVoxel,iiVoxel).sig] = getCausality(reshaped_MOT4,nTR,nTrials);

        disp('Causality - Low Attention');
        reshaped_MOT2(1,:,:) = area_MOT2_voxel(iVoxel,:,:);
        reshaped_MOT2(2,:,:) = area_MOT2_voxel(iiVoxel,:,:);
        
        [Granger_MOT2(iVoxel,iiVoxel).A,Granger_MOT2(iVoxel,iiVoxel).SIG,Granger_MOT2(iVoxel,iiVoxel).G,Granger_MOT2(iVoxel,iiVoxel).info,Granger_MOT2(iVoxel,iiVoxel).F,Granger_MOT2(iVoxel,iiVoxel).pval,Granger_MOT2(iVoxel,iiVoxel).sig] = getCausality(reshaped_MOT2,nTR,nTrials);

        disp('Causality - Resting State');
        reshaped_RestingState(1,:,:) = area_RestingState_voxel(iVoxel,:,:);
        reshaped_RestingState(2,:,:) = area_RestingState_voxel(iiVoxel,:,:);
        
        [Granger_RestingState(iVoxel,iiVoxel).A,Granger_RestingState(iVoxel,iiVoxel).SIG,Granger_RestingState(iVoxel,iiVoxel).G,Granger_RestingState(iVoxel,iiVoxel).info,Granger_RestingState(iVoxel,iiVoxel).F,Granger_RestingState(iVoxel,iiVoxel).pval,Granger_RestingState(iVoxel,iiVoxel).sig] = getCausality(reshaped_RestingState,nTR,nTrials);

    end
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','mvgc_AAL_DTI_3D','-',filename,'.mat'),'Granger_MOT4','Granger_MOT2','Granger_RestingState','-v7.3');

end

function [AIC,BIC,moAIC,moBIC,A,SIG,G,info,F,pval,sig] = getCausality(X,nTR,nTrials)

%% Parameters

ntrials   = nTrials;     % number of trials
nobs      = nTR;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

%morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
morder = 2;
momax     = 10;     % maximum model order for model order estimation

%acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
acmaxlags = 0;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

nvars = size(X,1); % number of variables

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% % Calculate information criteria up to specified maximum model order.
% 
% ptic('\n*** tsdata_to_infocrit\n');
% [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
% ptoc('*** tsdata_to_infocrit took ');
% 
% %amo = size(AT,3); % actual model order
% amo = 3;
% 
% fprintf('\nbest model order (AIC) = %d\n',moAIC);
% fprintf('best model order (BIC) = %d\n',moBIC);
% fprintf('actual model order     = %d\n',amo);
% 
% % Select model order.
% 
% if  strcmpi(morder,'actual')
%     morder = amo;
%     fprintf('\nusing actual model order = %d\n',morder);
% elseif strcmpi(morder,'AIC')
%     morder = moAIC;
%     fprintf('\nusing AIC best model order = %d\n',morder);
% elseif strcmpi(morder,'BIC')
%     morder = moBIC;
%     fprintf('\nusing BIC best model order = %d\n',morder);
% else
%     fprintf('\nusing specified model order = %d\n',morder);
% end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

end

function plotResults(settings,filename)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','mvgc_AAL_DTI_3D','-',filename,'.mat'));

%for iTrial=1:3
for iTrial=1:1
    
    %%% PLOT Causality Matrix

    alpha = 0.05;

    CLOW = min([min(min(Granger_MOT4(iTrial).F)), min(min(Granger_MOT2(iTrial).F)), min(min(Granger_RestingState(iTrial).F))]);
    CHIGH = max([max(max(Granger_MOT4(iTrial).F)), max(max(Granger_MOT2(iTrial).F)), max(max(Granger_RestingState(iTrial).F))]);

    CLIM = [CLOW CHIGH];

    f = figure;

    subplot(1,3,1);
    plot_pw(Granger_MOT4(iTrial).F,jet);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(Granger_MOT4(iTrial).pval,jet);
    title('p-values');
    subplot(1,3,3);
    plot_pw(Granger_MOT4(iTrial).sig,jet);
    title(['Significant at p = ' num2str(alpha)])

    %imagesc(rho_MOT4,CLIM);
    %colorbar;
    %title('Granger Causality');
    xlabel('AAL regions');
    ylabel('AAL regions');

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-HighAttention','.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-HighAttention','.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-HighAttention','.pdf'));

    g = figure;

    subplot(1,3,1);
    plot_pw(Granger_MOT2(iTrial).F,jet);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(Granger_MOT2(iTrial).pval,jet);
    title('p-values');
    subplot(1,3,3);
    plot_pw(Granger_MOT2(iTrial).sig,jet);
    title(['Significant at p = ' num2str(alpha)])

    %imagesc(rho_MOT2,CLIM);
    %colorbar;
    %title('Granger Causality');
    xlabel('AAL regions');
    ylabel('AAL regions');

    print(g,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-LowAttention','.jpg'));
    print(g,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-LowAttention','.eps'));
    print(g,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-LowAttention','.pdf'));

    h = figure;

    subplot(1,3,1);
    plot_pw(Granger_RestingState(iTrial).F,jet);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(Granger_RestingState(iTrial).pval,jet);
    title('p-values');
    subplot(1,3,3);
    plot_pw(Granger_RestingState(iTrial).sig,jet);
    title(['Significant at p = ' num2str(alpha)])

    %imagesc(rho_RestingState,CLIM);
    %colorbar;
    %title('Granger Causality');
    xlabel('AAL regions');
    ylabel('AAL regions');

    print(h,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-RestingState','.jpg'));
    print(h,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-RestingState','.eps'));
    print(h,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Granger-RestingState','.pdf'));

end

end

function reshaped = myReshapeRun(run,n)


nTrials = n;

nVoxels = size(run,1);
nTR = size(run,2);

reshaped = zeros(nVoxels,nTR/nTrials,nTrials);

for iTrial=1:nTrials
    
    for iVoxel=1:nVoxels
        
        reshaped(iVoxel,:,iTrial) = run(iVoxel,(1 + (nTR/nTrials)*(iTrial-1)):(nTR/nTrials + (nTR/nTrials)*(iTrial-1)));
        
    end
    
end

end

function processed_voxel  = preprocessVoxel(voxel)

    z_voxel = zscore(voxel);
    
    voxel_detrended = detrend(z_voxel);

    processed_voxel = voxel_detrended - mean(voxel_detrended);

end
