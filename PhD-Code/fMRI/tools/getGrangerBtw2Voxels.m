function [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(Voxel1,Voxel2)

nTR = length(Voxel1);
nvars = 2;

if size(Voxel1,2) == 1; Voxel1 = Voxel1'; end
if size(Voxel2,2) == 1; Voxel2 = Voxel2'; end

X = [Voxel1;Voxel2];

%% Parameters

% ntrials   = 10;     % number of trials
% nobs      = 1000;   % number of observations per trial

ntrials = 1;
nobs = nTR;

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
morder = 2;

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

A = nan(nvars,nvars);
SIG = nan(nvars,nvars); 
G = nan(nvars,nvars);
info = []; 
GF = nan(nvars,nvars);
Gpval = nan(nvars,nvars);
GSig = nan(nvars,nvars);

%% Model order estimation (<mvgc_schema.html#3 |A2|>) % Calculate information criteria up to specified maximum model order.

%[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);

%morder = moBIC;

%% VAR model estimation (<mvgc_schema.html#3 |A2|>) % Estimate VAR model of selected order from data.

[A,SIG] = tsdata_to_var(X,morder,regmode);

% Check for failed regression

%assert(~isbad(A),'VAR estimation failed');

if ~isbad(A)
    
    %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

    [G,info] = var_to_autocov(A,SIG,acmaxlags);

    %var_info(info,true); % report results (and bail out on error)

    if ~info.error
        
        %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

        % Calculate time-domain pairwise-conditional causalities - this just requires
        % the autocovariance sequence.

        GF = autocov_to_pwcgc(G);

        if ~isbad(GF,false)

            %assert(~isbad(GF,false),'GC calculation failed');

            % Significance test using theoretical null distribution, adjusting for multiple
            % hypotheses.

            Gpval = mvgc_pval(GF,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
            GSig  = significance(Gpval,alpha,mhtc);

        end
    
    end
    
end

end

