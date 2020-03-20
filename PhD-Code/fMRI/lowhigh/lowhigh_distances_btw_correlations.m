function lowhigh_distances_btw_correlations

settings_jan_0805;
%settings_elena_2905;

%doTheMath(settings);

% doTheMathMean(settings);
% doTheMathPICA(settings);
% doTheMathVoxel(settings);

%plotStuff(settings);

%plotStuffbyRun(settings,'Mean');
%plotStuffbyRun(settings,'PICA');
plotStuffbyRun(settings,'Voxel');

end

function doTheMath(settings)

step = 0.0001;
interval = -1:step:1;
nDistributions = 24;

all_distributions = zeros(nDistributions,length(interval));

disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM MEAN OF VOXELS CORRELATIONS');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'));

rho_MOT4Run1_FPO = getMeanRhoFPO(rho_MOT4Run1);
rho_MOT4Run2_FPO = getMeanRhoFPO(rho_MOT4Run2);

rho_RestingStateRun1_FPO = getMeanRhoFPO(rho_RestingStateRun1);
rho_RestingStateRun2_FPO = getMeanRhoFPO(rho_RestingStateRun2);

pval_MOT4Run1_FPO = getMeanRhoFPO(pval_MOT4Run1);
pval_MOT4Run2_FPO = getMeanRhoFPO(pval_MOT4Run2);

pval_RestingStateRun1_FPO = getMeanRhoFPO(pval_RestingStateRun1);
pval_RestingStateRun2_FPO = getMeanRhoFPO(pval_RestingStateRun2);

all_distributions(1,:) = getProbabilityDistribution(rho_MOT4Run1_FPO,pval_MOT4Run1_FPO);
all_distributions(2,:) = getProbabilityDistribution(rho_MOT4Run2_FPO,pval_MOT4Run2_FPO);

distributions_labels{1} = 'High-Run-1-Mean';
distributions_labels{2} = 'High-Run-2-Mean';

all_distributions(19,:) = getProbabilityDistribution(rho_RestingStateRun1_FPO,pval_RestingStateRun1_FPO);
all_distributions(20,:) = getProbabilityDistribution(rho_RestingStateRun2_FPO,pval_RestingStateRun2_FPO);

distributions_labels{19} = 'RestingState-Run-1-Mean';
distributions_labels{20} = 'RestingState-Run-2-Mean';

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO-Half','.mat'));

% rho_MOT4Run1_FPO_FH = getMeanRhoFPO(rho_MOT4Run1_FH);
% rho_MOT4Run2_FPO_FH = getMeanRhoFPO(rho_MOT4Run2_FH);
% rho_MOT4Run1_FPO_SH = getMeanRhoFPO(rho_MOT4Run1_SH);
% rho_MOT4Run2_FPO_SH = getMeanRhoFPO(rho_MOT4Run2_SH);
% 
% rho_RestingStateRun1_FPO_FH = getMeanRhoFPO(rho_RestingStateRun1_FH);
% rho_RestingStateRun2_FPO_FH = getMeanRhoFPO(rho_RestingStateRun2_FH);
% rho_RestingStateRun1_FPO_SH = getMeanRhoFPO(rho_RestingStateRun1_SH);
% rho_RestingStateRun2_FPO_SH = getMeanRhoFPO(rho_RestingStateRun2_SH);
% 
% pval_MOT4Run1_FPO_FH = getMeanRhoFPO(pval_MOT4Run1_FH);
% pval_MOT4Run2_FPO_FH = getMeanRhoFPO(pval_MOT4Run2_FH);
% pval_MOT4Run1_FPO_SH = getMeanRhoFPO(pval_MOT4Run1_SH);
% pval_MOT4Run2_FPO_SH = getMeanRhoFPO(pval_MOT4Run2_SH);
% 
% pval_RestingStateRun1_FPO_FH = getMeanRhoFPO(pval_RestingStateRun1_FH);
% pval_RestingStateRun2_FPO_FH = getMeanRhoFPO(pval_RestingStateRun2_FH);
% pval_RestingStateRun1_FPO_SH = getMeanRhoFPO(pval_RestingStateRun1_SH);
% pval_RestingStateRun2_FPO_SH = getMeanRhoFPO(pval_RestingStateRun2_SH);

all_distributions(3,:) = getProbabilityDistribution(rho_MOT4Run1_FH,pval_MOT4Run1_FH);
all_distributions(4,:) = getProbabilityDistribution(rho_MOT4Run2_FH,pval_MOT4Run2_FH);
all_distributions(5,:) = getProbabilityDistribution(rho_MOT4Run1_SH,pval_MOT4Run1_SH);
all_distributions(6,:) = getProbabilityDistribution(rho_MOT4Run2_SH,pval_MOT4Run2_SH);

distributions_labels{3} = 'High-Run-1-Mean-FH';
distributions_labels{4} = 'High-Run-2-Mean-FH';
distributions_labels{5} = 'High-Run-1-Mean-SH';
distributions_labels{6} = 'High-Run-2-Mean-SH';

all_distributions(21,:) = getProbabilityDistribution(rho_RestingStateRun1_FH,pval_RestingStateRun1_FH);
all_distributions(22,:) = getProbabilityDistribution(rho_RestingStateRun2_FH,pval_RestingStateRun2_FH);
all_distributions(23,:) = getProbabilityDistribution(rho_RestingStateRun1_SH,pval_RestingStateRun1_SH);
all_distributions(24,:) = getProbabilityDistribution(rho_RestingStateRun2_SH,pval_RestingStateRun2_SH);

distributions_labels{21} = 'Rest-Run-1-Mean-FH';
distributions_labels{22} = 'Rest-Run-2-Mean-FH';
distributions_labels{23} = 'Rest-Run-1-Mean-SH';
distributions_labels{24} = 'Rest-Run-2-Mean-SH';

clear rho_MOT4Run1_FH
clear rho_MOT4Run2_FH
clear rho_MOT4Run1_SH
clear rho_MOT4Run2_SH
clear rho_RestingStateRun1_FH
clear rho_RestingStateRun2_FH
clear rho_RestingStateRun1_SH
clear rho_RestingStateRun2_SH

clear rho_MOT4Run1_FPO_FH
clear rho_MOT4Run2_FPO_FH
clear rho_MOT4Run1_FPO_SH
clear rho_MOT4Run2_FPO_SH
clear rho_RestingStateRun1_FPO_FH
clear rho_RestingStateRun2_FPO_FH
clear rho_RestingStateRun1_FPO_SH
clear rho_RestingStateRun2_FPO_SH


disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM PICA CORRELATIONS');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

FF = FC_ICA_FO.run(1).rho_high(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FO','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FP','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-PO','.mat'));

FF = FC_ICA_FO.run(1).rho_high_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_high_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

all_distributions(7,:) = getProbabilityDistribution(rho_high_run1_FPO,pval_high_run1_FPO);
all_distributions(8,:) = getProbabilityDistribution(rho_high_run2_FPO,pval_high_run2_FPO);

distributions_labels{7} = 'High-Run-1-PICA';
distributions_labels{8} = 'High-Run-2-PICA';

all_distributions(25,:) = getProbabilityDistribution(rho_rest_run1_FPO,pval_rest_run1_FPO);
all_distributions(26,:) = getProbabilityDistribution(rho_rest_run2_FPO,pval_rest_run2_FPO);

distributions_labels{25} = 'Rest-Run-1-PICA';
distributions_labels{26} = 'Rest-Run-2-PICA';

all_distributions(9,:) = getProbabilityDistribution(rho_high_run1_FPO_FH,pval_high_run1_FPO_FH);
all_distributions(10,:) = getProbabilityDistribution(rho_high_run2_FPO_FH,pval_high_run2_FPO_FH);
all_distributions(11,:) = getProbabilityDistribution(rho_high_run1_FPO_SH,pval_high_run1_FPO_SH);
all_distributions(12,:) = getProbabilityDistribution(rho_high_run2_FPO_SH,pval_high_run2_FPO_SH);

distributions_labels{9} = 'High-Run-1-PICA-FH';
distributions_labels{10} = 'High-Run-2-PICA-FH';
distributions_labels{11} = 'High-Run-1-PICA-SH';
distributions_labels{12} = 'High-Run-2-PICA-SH';

all_distributions(27,:) = getProbabilityDistribution(rho_rest_run1_FPO_FH,pval_rest_run1_FPO_FH);
all_distributions(28,:) = getProbabilityDistribution(rho_rest_run2_FPO_FH,pval_rest_run2_FPO_FH);
all_distributions(29,:) = getProbabilityDistribution(rho_rest_run1_FPO_SH,pval_rest_run1_FPO_SH);
all_distributions(30,:) = getProbabilityDistribution(rho_rest_run2_FPO_SH,pval_rest_run2_FPO_SH);

distributions_labels{27} = 'Rest-Run-1-PICA-FH';
distributions_labels{28} = 'Rest-Run-2-PICA-FH';
distributions_labels{29} = 'Rest-Run-1-PICA-SH';
distributions_labels{30} = 'Rest-Run-2-PICA-SH';

clear rho_high_run1_FPO_FH
clear rho_high_run2_FPO_FH
clear rho_high_run1_FPO_SH
clear rho_high_run2_FPO_SH
clear rho_rest_run1_FPO_FH
clear rho_rest_run2_FPO_FH
clear rho_rest_run1_FPO_SH
clear rho_rest_run2_FPO_SH

disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM VOXEL LEVEL CORRELATIONS');

load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;
idx_ROI = 66;

label_ROI = AAL_ROI(idx_ROI).Nom_L;
label_ROI = strrep(label_ROI,'_','-');
disp(label_ROI);

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',label_ROI,'.mat'));

rho_high_run1_voxel = FC_Voxels.run(1).rho_MOT4;
rho_high_run2_voxel = FC_Voxels.run(2).rho_MOT4;

rho_rest_run1_voxel = FC_Voxels.run(1).rho_RestingState;
rho_rest_run2_voxel = FC_Voxels.run(2).rho_RestingState;

pval_high_run1_voxel = FC_Voxels.run(1).pval_MOT4;
pval_high_run2_voxel = FC_Voxels.run(2).pval_MOT4;

pval_rest_run1_voxel = FC_Voxels.run(1).pval_RestingState;
pval_rest_run2_voxel = FC_Voxels.run(2).pval_RestingState;

all_distributions(13,:) = getProbabilityDistribution(rho_high_run1_voxel,pval_high_run1_voxel);
all_distributions(14,:) = getProbabilityDistribution(rho_high_run2_voxel,pval_high_run2_voxel);

distributions_labels{13} = 'High-Run-1-Voxel-FH';
distributions_labels{14} = 'High-Run-2-Voxel-FH';

all_distributions(31,:) = getProbabilityDistribution(rho_rest_run1_voxel,pval_rest_run1_voxel);
all_distributions(32,:) = getProbabilityDistribution(rho_rest_run2_voxel,pval_rest_run2_voxel);

distributions_labels{31} = 'Rest-Run-1-Voxel-FH';
distributions_labels{32} = 'Rest-Run-2-Voxel-FH';

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',label_ROI,'.mat'));

rho_high_run1_voxel_FH = FC_Voxels.run(1).rho_MOT4_FH;
rho_high_run2_voxel_FH = FC_Voxels.run(2).rho_MOT4_FH;
rho_high_run1_voxel_SH = FC_Voxels.run(1).rho_MOT4_SH;
rho_high_run2_voxel_SH = FC_Voxels.run(2).rho_MOT4_SH;
rho_rest_run1_voxel_FH = FC_Voxels.run(1).rho_RestingState_FH;
rho_rest_run2_voxel_FH = FC_Voxels.run(2).rho_RestingState_FH;
rho_rest_run1_voxel_SH = FC_Voxels.run(1).rho_RestingState_SH;
rho_rest_run2_voxel_SH = FC_Voxels.run(2).rho_RestingState_SH;

pval_high_run1_voxel_FH = FC_Voxels.run(1).pval_MOT4_FH;
pval_high_run2_voxel_FH = FC_Voxels.run(2).pval_MOT4_FH;
pval_high_run1_voxel_SH = FC_Voxels.run(1).pval_MOT4_SH;
pval_high_run2_voxel_SH = FC_Voxels.run(2).pval_MOT4_SH;
pval_rest_run1_voxel_FH = FC_Voxels.run(1).pval_RestingState_FH;
pval_rest_run2_voxel_FH = FC_Voxels.run(2).pval_RestingState_FH;
pval_rest_run1_voxel_SH = FC_Voxels.run(1).pval_RestingState_SH;
pval_rest_run2_voxel_SH = FC_Voxels.run(2).pval_RestingState_SH;

all_distributions(15,:) = getProbabilityDistribution(rho_high_run1_voxel_FH,pval_high_run1_voxel_FH);
all_distributions(16,:) = getProbabilityDistribution(rho_high_run2_voxel_FH,pval_high_run2_voxel_FH);
all_distributions(17,:) = getProbabilityDistribution(rho_high_run1_voxel_SH,pval_high_run1_voxel_SH);
all_distributions(18,:) = getProbabilityDistribution(rho_high_run2_voxel_SH,pval_high_run2_voxel_SH);

distributions_labels{15} = 'High-Run-1-Voxel-FH';
distributions_labels{16} = 'High-Run-2-Voxel-FH';
distributions_labels{17} = 'High-Run-1-Voxel-SH';
distributions_labels{18} = 'High-Run-2-Voxel-SH';

all_distributions(33,:) = getProbabilityDistribution(rho_rest_run1_voxel_FH,pval_rest_run1_voxel_FH);
all_distributions(34,:) = getProbabilityDistribution(rho_rest_run2_voxel_FH,pval_rest_run2_voxel_FH);
all_distributions(35,:) = getProbabilityDistribution(rho_rest_run1_voxel_SH,pval_rest_run1_voxel_SH);
all_distributions(36,:) = getProbabilityDistribution(rho_rest_run2_voxel_SH,pval_rest_run2_voxel_SH);

distributions_labels{33} = 'Rest-Run-1-Voxel-FH';
distributions_labels{34} = 'Rest-Run-2-Voxel-FH';
distributions_labels{35} = 'Rest-Run-1-Voxel-SH';
distributions_labels{36} = 'Rest-Run-2-Voxel-SH';

clear rho_high_run1_voxel_FH
clear rho_high_run2_voxel_FH
clear rho_high_run1_voxel_SH
clear rho_high_run2_voxel_SH
clear rho_rest_run1_voxel_FH
clear rho_rest_run2_voxel_FH
clear rho_rest_run1_voxel_SH
clear rho_rest_run2_voxel_SH

distances = getDistances(all_distributions);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances','.mat'),'all_distributions','distributions_labels','distances');

end

function doTheMathMean(settings)

step = 0.001;
interval = -1:step:1;
nDistributions = 12;

all_distributions = zeros(nDistributions,length(interval));

disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM MEAN OF VOXELS CORRELATIONS');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'));

rho_MOT4Run1_FPO = getMeanRhoFPO(rho_MOT4Run1);
rho_MOT4Run2_FPO = getMeanRhoFPO(rho_MOT4Run2);

rho_RestingStateRun1_FPO = getMeanRhoFPO(rho_RestingStateRun1);
rho_RestingStateRun2_FPO = getMeanRhoFPO(rho_RestingStateRun2);

pval_MOT4Run1_FPO = getMeanRhoFPO(pval_MOT4Run1);
pval_MOT4Run2_FPO = getMeanRhoFPO(pval_MOT4Run2);

pval_RestingStateRun1_FPO = getMeanRhoFPO(pval_RestingStateRun1);
pval_RestingStateRun2_FPO = getMeanRhoFPO(pval_RestingStateRun2);

all_distributions(1,:) = getProbabilityDistribution(rho_MOT4Run1_FPO,pval_MOT4Run1_FPO);
all_distributions(2,:) = getProbabilityDistribution(rho_MOT4Run2_FPO,pval_MOT4Run2_FPO);

distributions_labels{1} = 'High-Run-1-Mean';
distributions_labels{2} = 'High-Run-2-Mean';

all_distributions(7,:) = getProbabilityDistribution(rho_RestingStateRun1_FPO,pval_RestingStateRun1_FPO);
all_distributions(8,:) = getProbabilityDistribution(rho_RestingStateRun2_FPO,pval_RestingStateRun2_FPO);

distributions_labels{7} = 'RestingState-Run-1-Mean';
distributions_labels{8} = 'RestingState-Run-2-Mean';

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO-Half','.mat'));

% rho_MOT4Run1_FPO_FH = getMeanRhoFPO(rho_MOT4Run1_FH);
% rho_MOT4Run2_FPO_FH = getMeanRhoFPO(rho_MOT4Run2_FH);
% rho_MOT4Run1_FPO_SH = getMeanRhoFPO(rho_MOT4Run1_SH);
% rho_MOT4Run2_FPO_SH = getMeanRhoFPO(rho_MOT4Run2_SH);
% 
% rho_RestingStateRun1_FPO_FH = getMeanRhoFPO(rho_RestingStateRun1_FH);
% rho_RestingStateRun2_FPO_FH = getMeanRhoFPO(rho_RestingStateRun2_FH);
% rho_RestingStateRun1_FPO_SH = getMeanRhoFPO(rho_RestingStateRun1_SH);
% rho_RestingStateRun2_FPO_SH = getMeanRhoFPO(rho_RestingStateRun2_SH);
% 
% pval_MOT4Run1_FPO_FH = getMeanRhoFPO(pval_MOT4Run1_FH);
% pval_MOT4Run2_FPO_FH = getMeanRhoFPO(pval_MOT4Run2_FH);
% pval_MOT4Run1_FPO_SH = getMeanRhoFPO(pval_MOT4Run1_SH);
% pval_MOT4Run2_FPO_SH = getMeanRhoFPO(pval_MOT4Run2_SH);
% 
% pval_RestingStateRun1_FPO_FH = getMeanRhoFPO(pval_RestingStateRun1_FH);
% pval_RestingStateRun2_FPO_FH = getMeanRhoFPO(pval_RestingStateRun2_FH);
% pval_RestingStateRun1_FPO_SH = getMeanRhoFPO(pval_RestingStateRun1_SH);
% pval_RestingStateRun2_FPO_SH = getMeanRhoFPO(pval_RestingStateRun2_SH);

all_distributions(3,:) = getProbabilityDistribution(rho_MOT4Run1_FH,pval_MOT4Run1_FH);
all_distributions(4,:) = getProbabilityDistribution(rho_MOT4Run2_FH,pval_MOT4Run2_FH);
all_distributions(5,:) = getProbabilityDistribution(rho_MOT4Run1_SH,pval_MOT4Run1_SH);
all_distributions(6,:) = getProbabilityDistribution(rho_MOT4Run2_SH,pval_MOT4Run2_SH);

distributions_labels{3} = 'High-Run-1-Mean-FH';
distributions_labels{4} = 'High-Run-2-Mean-FH';
distributions_labels{5} = 'High-Run-1-Mean-SH';
distributions_labels{6} = 'High-Run-2-Mean-SH';

all_distributions(9,:) = getProbabilityDistribution(rho_RestingStateRun1_FH,pval_RestingStateRun1_FH);
all_distributions(10,:) = getProbabilityDistribution(rho_RestingStateRun2_FH,pval_RestingStateRun2_FH);
all_distributions(11,:) = getProbabilityDistribution(rho_RestingStateRun1_SH,pval_RestingStateRun1_SH);
all_distributions(12,:) = getProbabilityDistribution(rho_RestingStateRun2_SH,pval_RestingStateRun2_SH);

distributions_labels{9} = 'Rest-Run-1-Mean-FH';
distributions_labels{10} = 'Rest-Run-2-Mean-FH';
distributions_labels{11} = 'Rest-Run-1-Mean-SH';
distributions_labels{12} = 'Rest-Run-2-Mean-SH';

clear rho_MOT4Run1_FH
clear rho_MOT4Run2_FH
clear rho_MOT4Run1_SH
clear rho_MOT4Run2_SH
clear rho_RestingStateRun1_FH
clear rho_RestingStateRun2_FH
clear rho_RestingStateRun1_SH
clear rho_RestingStateRun2_SH

clear rho_MOT4Run1_FPO_FH
clear rho_MOT4Run2_FPO_FH
clear rho_MOT4Run1_FPO_SH
clear rho_MOT4Run2_FPO_SH
clear rho_RestingStateRun1_FPO_FH
clear rho_RestingStateRun2_FPO_FH
clear rho_RestingStateRun1_FPO_SH
clear rho_RestingStateRun2_FPO_SH

distances = getDistances(all_distributions);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-Mean','.mat'),'all_distributions','distributions_labels','distances');

end

function doTheMathPICA(settings)

step = 0.001;
interval = -1:step:1;
nDistributions = 12;

all_distributions = zeros(nDistributions,length(interval));

disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM PICA CORRELATIONS');

disp('loading full');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe','.mat'));

FF = FC_ICA_FO.run(1).rho_high(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO = [FF(:);FO(:);FP(:);PO(:)];

disp('loading half');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FO','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-FP','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','groupICA-ROI_grouped_by_lobe_Half-PO','.mat'));

FF = FC_ICA_FO.run(1).rho_high_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_high_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_high_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_high_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_high_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_high_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_high_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_high_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_high_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_high_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_high_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_high_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_high_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_high_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_high_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_high_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_high_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_high_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_high_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_high_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_high_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest_FH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest_FH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest_FH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest_FH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest_FH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest_FH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO_FH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).rho_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).rho_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).rho_rest_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).rho_rest_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
rho_rest_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(1).pval_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,1:FC_n_FO.run(1).nIC_area1_total);
FO = FC_ICA_FO.run(1).pval_rest_SH(1:FC_n_FO.run(1).nIC_area1_total,(FC_n_FO.run(1).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(1).pval_rest_SH(1:FC_n_FP.run(1).nIC_area1_total,(FC_n_FP.run(1).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(1).pval_rest_SH(1:FC_n_PO.run(1).nIC_area1_total,(FC_n_PO.run(1).nIC_area1_total+1):end);
pval_rest_run1_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).rho_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).rho_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).rho_rest_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).rho_rest_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
rho_rest_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

FF = FC_ICA_FO.run(2).pval_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,1:FC_n_FO.run(2).nIC_area1_total);
FO = FC_ICA_FO.run(2).pval_rest_SH(1:FC_n_FO.run(2).nIC_area1_total,(FC_n_FO.run(2).nIC_area1_total+1):end);
FP = FC_ICA_FP.run(2).pval_rest_SH(1:FC_n_FP.run(2).nIC_area1_total,(FC_n_FP.run(2).nIC_area1_total+1):end);
PO = FC_ICA_PO.run(2).pval_rest_SH(1:FC_n_PO.run(2).nIC_area1_total,(FC_n_PO.run(2).nIC_area1_total+1):end);
pval_rest_run2_FPO_SH = [FF(:);FO(:);FP(:);PO(:)];

disp('distribution 1,2');

all_distributions(1,:) = getProbabilityDistribution(rho_high_run1_FPO,pval_high_run1_FPO);
all_distributions(2,:) = getProbabilityDistribution(rho_high_run2_FPO,pval_high_run2_FPO);

distributions_labels{1} = 'High-Run-1-PICA';
distributions_labels{2} = 'High-Run-2-PICA';

disp('distribution 7,8');

all_distributions(7,:) = getProbabilityDistribution(rho_rest_run1_FPO,pval_rest_run1_FPO);
all_distributions(8,:) = getProbabilityDistribution(rho_rest_run2_FPO,pval_rest_run2_FPO);

distributions_labels{7} = 'Rest-Run-1-PICA';
distributions_labels{8} = 'Rest-Run-2-PICA';

disp('distribution 3,4,5,6');

all_distributions(3,:) = getProbabilityDistribution(rho_high_run1_FPO_FH,pval_high_run1_FPO_FH);
all_distributions(4,:) = getProbabilityDistribution(rho_high_run2_FPO_FH,pval_high_run2_FPO_FH);
all_distributions(5,:) = getProbabilityDistribution(rho_high_run1_FPO_SH,pval_high_run1_FPO_SH);
all_distributions(6,:) = getProbabilityDistribution(rho_high_run2_FPO_SH,pval_high_run2_FPO_SH);

distributions_labels{3} = 'High-Run-1-PICA-FH';
distributions_labels{4} = 'High-Run-2-PICA-FH';
distributions_labels{5} = 'High-Run-1-PICA-SH';
distributions_labels{6} = 'High-Run-2-PICA-SH';

disp('distribution 9,10,11,12');

all_distributions(9,:) = getProbabilityDistribution(rho_rest_run1_FPO_FH,pval_rest_run1_FPO_FH);
all_distributions(10,:) = getProbabilityDistribution(rho_rest_run2_FPO_FH,pval_rest_run2_FPO_FH);
all_distributions(11,:) = getProbabilityDistribution(rho_rest_run1_FPO_SH,pval_rest_run1_FPO_SH);
all_distributions(12,:) = getProbabilityDistribution(rho_rest_run2_FPO_SH,pval_rest_run2_FPO_SH);

distributions_labels{9} = 'Rest-Run-1-PICA-FH';
distributions_labels{10} = 'Rest-Run-2-PICA-FH';
distributions_labels{11} = 'Rest-Run-1-PICA-SH';
distributions_labels{12} = 'Rest-Run-2-PICA-SH';

clear rho_high_run1_FPO_FH
clear rho_high_run2_FPO_FH
clear rho_high_run1_FPO_SH
clear rho_high_run2_FPO_SH
clear rho_rest_run1_FPO_FH
clear rho_rest_run2_FPO_FH
clear rho_rest_run1_FPO_SH
clear rho_rest_run2_FPO_SH

distances = getDistances(all_distributions);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-PICA','.mat'),'all_distributions','distributions_labels','distances');

end

function doTheMathVoxel(settings)

step = 0.001;
interval = -1:step:1;
nDistributions = 8;

all_distributions = zeros(nDistributions,length(interval));

disp('GET FULL/HALF TIME SERIES DISTRIBUTIONS FROM VOXEL LEVEL CORRELATIONS');

load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;
idx_ROI = 66;

label_ROI = AAL_ROI(idx_ROI).Nom_L;
label_ROI = strrep(label_ROI,'_','-');
disp(label_ROI);

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',label_ROI,'.mat'));

rho_high_run1_voxel = FC_Voxels.run(1).rho_MOT4;
rho_high_run2_voxel = FC_Voxels.run(2).rho_MOT4;

rho_rest_run1_voxel = FC_Voxels.run(1).rho_RestingState;
rho_rest_run2_voxel = FC_Voxels.run(2).rho_RestingState;

pval_high_run1_voxel = FC_Voxels.run(1).pval_MOT4;
pval_high_run2_voxel = FC_Voxels.run(2).pval_MOT4;

pval_rest_run1_voxel = FC_Voxels.run(1).pval_RestingState;
pval_rest_run2_voxel = FC_Voxels.run(2).pval_RestingState;

all_distributions(1,:) = getProbabilityDistribution(rho_high_run1_voxel,pval_high_run1_voxel);
all_distributions(2,:) = getProbabilityDistribution(rho_high_run2_voxel,pval_high_run2_voxel);

distributions_labels{1} = 'High-Run-1-Voxel-FH';
distributions_labels{2} = 'High-Run-2-Voxel-FH';

all_distributions(7,:) = getProbabilityDistribution(rho_rest_run1_voxel,pval_rest_run1_voxel);
all_distributions(8,:) = getProbabilityDistribution(rho_rest_run2_voxel,pval_rest_run2_voxel);

distributions_labels{7} = 'Rest-Run-1-Voxel-FH';
distributions_labels{8} = 'Rest-Run-2-Voxel-FH';

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',label_ROI,'.mat'));

rho_high_run1_voxel_FH = FC_Voxels.run(1).rho_MOT4_FH;
rho_high_run2_voxel_FH = FC_Voxels.run(2).rho_MOT4_FH;
rho_high_run1_voxel_SH = FC_Voxels.run(1).rho_MOT4_SH;
rho_high_run2_voxel_SH = FC_Voxels.run(2).rho_MOT4_SH;
rho_rest_run1_voxel_FH = FC_Voxels.run(1).rho_RestingState_FH;
rho_rest_run2_voxel_FH = FC_Voxels.run(2).rho_RestingState_FH;
rho_rest_run1_voxel_SH = FC_Voxels.run(1).rho_RestingState_SH;
rho_rest_run2_voxel_SH = FC_Voxels.run(2).rho_RestingState_SH;

pval_high_run1_voxel_FH = FC_Voxels.run(1).pval_MOT4_FH;
pval_high_run2_voxel_FH = FC_Voxels.run(2).pval_MOT4_FH;
pval_high_run1_voxel_SH = FC_Voxels.run(1).pval_MOT4_SH;
pval_high_run2_voxel_SH = FC_Voxels.run(2).pval_MOT4_SH;
pval_rest_run1_voxel_FH = FC_Voxels.run(1).pval_RestingState_FH;
pval_rest_run2_voxel_FH = FC_Voxels.run(2).pval_RestingState_FH;
pval_rest_run1_voxel_SH = FC_Voxels.run(1).pval_RestingState_SH;
pval_rest_run2_voxel_SH = FC_Voxels.run(2).pval_RestingState_SH;

all_distributions(3,:) = getProbabilityDistribution(rho_high_run1_voxel_FH,pval_high_run1_voxel_FH);
all_distributions(4,:) = getProbabilityDistribution(rho_high_run2_voxel_FH,pval_high_run2_voxel_FH);
all_distributions(5,:) = getProbabilityDistribution(rho_high_run1_voxel_SH,pval_high_run1_voxel_SH);
all_distributions(6,:) = getProbabilityDistribution(rho_high_run2_voxel_SH,pval_high_run2_voxel_SH);

distributions_labels{3} = 'High-Run-1-Voxel-FH';
distributions_labels{4} = 'High-Run-2-Voxel-FH';
distributions_labels{5} = 'High-Run-1-Voxel-SH';
distributions_labels{6} = 'High-Run-2-Voxel-SH';

all_distributions(9,:) = getProbabilityDistribution(rho_rest_run1_voxel_FH,pval_rest_run1_voxel_FH);
all_distributions(10,:) = getProbabilityDistribution(rho_rest_run2_voxel_FH,pval_rest_run2_voxel_FH);
all_distributions(11,:) = getProbabilityDistribution(rho_rest_run1_voxel_SH,pval_rest_run1_voxel_SH);
all_distributions(12,:) = getProbabilityDistribution(rho_rest_run2_voxel_SH,pval_rest_run2_voxel_SH);

distributions_labels{9} = 'Rest-Run-1-Voxel-FH';
distributions_labels{10} = 'Rest-Run-2-Voxel-FH';
distributions_labels{11} = 'Rest-Run-1-Voxel-SH';
distributions_labels{12} = 'Rest-Run-2-Voxel-SH';

clear rho_high_run1_voxel_FH
clear rho_high_run2_voxel_FH
clear rho_high_run1_voxel_SH
clear rho_high_run2_voxel_SH
clear rho_rest_run1_voxel_FH
clear rho_rest_run2_voxel_FH
clear rho_rest_run1_voxel_SH
clear rho_rest_run2_voxel_SH

distances = getDistances(all_distributions);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-Voxel-',label_ROI,'.mat'),'all_distributions','distributions_labels','distances');

end

function plotStuff(settings)

%load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances','.mat'));
load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-0.001-Angular-R','.mat'));

label = 'High-Mean';
dist = all_distributions(1:4,:);
plot4Distributions(settings,dist,label);

label = 'High-PICA';
dist = all_distributions(5:8,:);
plot4Distributions(settings,dist,label);

label = 'High-Angular-R-Voxels';
dist = all_distributions(9:12,:);
plot4Distributions(settings,dist,label);

label = 'Rest-Mean';
dist = all_distributions(13:16,:);
plot4Distributions(settings,dist,label);

label = 'Rest-PICA';
dist = all_distributions(17:20,:);
plot4Distributions(settings,dist,label);

label = 'Rest-Angular-R-Voxels';
dist = all_distributions(21:24,:);
plot4Distributions(settings,dist,label);

end

function probability_distribution = getProbabilityDistribution(rho_mat,pval_mat)

pcriterion = 0.01;
idx = find(pval_mat>pcriterion);
rho_mat(idx) = 0;

rho_mat_nozeros = rho_mat;
rho_mat_nozeros(rho_mat==0) = [];

step = 0.001;
interval = -1:step:1;

rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

k = ksdensity(rho_mat_fisher(:),interval);

probability_distribution = k .* step;

end

function mean_rho_FPO = getMeanRhoFPO(rho_mat)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
aal_idx = [idx_frontal,idx_occipital,idx_parietal];

new_rho = zeros(size(aal_idx,2),size(aal_idx,2));

for iidx=1:length(aal_idx)
    
    new_rho(iidx,:) = rho_mat(aal_idx(iidx),aal_idx);
    
end

mean_rho_FPO = new_rho;

end

function distance = getDistances(all_distributions)

nDistributions = size(all_distributions,1);

distance = zeros(nDistributions);

for iDis=1:nDistributions
   
    X = all_distributions(iDis,:);
    
    for iiDis=iDis:nDistributions
        
        Y = all_distributions(iiDis,:);
        
        XY = X.*Y;
        sqrtXY = sqrt(XY);
        
        distance(iDis,iiDis) = sqrt( 1 - sum( sqrtXY ) );
        distance(iiDis,iDis) = distance(iDis,iiDis);
        
    end
    
end

end

function plot4Distributions(settings,dist,label)

f = figure;

interval = -1:0.001:1;

plot(interval,dist(1,:),'b');
hold on
plot(interval,dist(2,:),'r');
plot(interval,dist(3,:),'y');
plot(interval,dist(4,:),'g');

ylim([0 0.005]);
xlim([-1 1]);
title(label);
legend({'Run1-FH' 'Run2-FH' 'Run1-SH' 'Run2-SH'});
xlabel('Pearson Correlation');
ylabel('Probability');

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.jpg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.eps'));
print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.pdf'));

end

function plotStuffbyRun(settings,kind)

if strcmp(kind,'Voxel')
    
    load_roi = load('ROI_MNI_V4_List.mat');
    AAL_ROI = load_roi.ROI;
    idx_ROI = 66;

    label_ROI = AAL_ROI(idx_ROI).Nom_L;
    label_ROI = strrep(label_ROI,'_','-');
    disp(label_ROI);
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-Voxel-',label_ROI,'.mat'));

else
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','Distributions-Distances-',kind,'.mat'));

end


label = strcat('High-',kind);
dist = all_distributions(1:6,:);
plot6Distributions(settings,dist,label);

label = strcat('Rest-',kind);
dist = all_distributions(7:12,:);
plot6Distributions(settings,dist,label);


end

function plot6Distributions(settings,dist,label)

f = figure;

interval = -1:0.001:1;

plot(interval,dist(1,:),'b');
hold on
plot(interval,dist(2,:),'r');
plot(interval,dist(3,:),'y');
plot(interval,dist(4,:),'g');
plot(interval,dist(5,:),'k');
plot(interval,dist(6,:),'m');

ylim([0 0.005]);
xlim([-1 1]);
title(label);
legend({'Run1-Full' 'Run2-Full' 'Run1-FirstHalf' 'Run2-FirstHalf' 'Run1-SecondHalf' 'Run2-SecondHalf'});
xlabel('Pearson Correlation');
ylabel('Probability');

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.jpg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.eps'));
print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Distributions-Distances-',label,'.pdf'));

end


