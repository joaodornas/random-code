
function plot_subj2_granger

Condition = 'RestingState';

Net1 = 'DAN';
Net2 = 'VAN';
[results_RestingState_DAN, results_RestingState_VAN, results_RestingState_DAN_VAN] = getResults(Condition,Net1,Net2);

Net1 = 'VAN';
Net2 = 'VIS';
[results_RestingState_VAN, results_RestingState_VIS, results_RestingState_VAN_VIS] = getResults(Condition,Net1,Net2);

Net1 = 'DAN';
Net2 = 'VIS';
[results_RestingState_DAN, results_RestingState_VIS, results_RestingState_DAN_VIS] = getResults(Condition,Net1,Net2);

plotCondition(Condition,results_RestingState_DAN,results_RestingState_VAN,results_RestingState_VIS,results_RestingState_DAN_VAN,results_RestingState_VAN_VIS,results_RestingState_DAN_VIS);

Condition = 'Passive';

Net1 = 'DAN';
Net2 = 'VAN';
[results_Passive_DAN, results_Passive_VAN, results_Passive_DAN_VAN] = getResults(Condition,Net1,Net2);

Net1 = 'VAN';
Net2 = 'VIS';
[results_Passive_VAN, results_Passive_VIS, results_Passive_VAN_VIS] = getResults(Condition,Net1,Net2);

Net1 = 'DAN';
Net2 = 'VIS';
[results_Passive_DAN, results_Passive_VIS, results_Passive_DAN_VIS] = getResults(Condition,Net1,Net2);

plotCondition(Condition,results_Passive_DAN,results_Passive_VAN,results_Passive_VIS,results_Passive_DAN_VAN,results_Passive_VAN_VIS,results_Passive_DAN_VIS);

Condition = 'Track';

Net1 = 'DAN';
Net2 = 'VAN';
[results_Track_DAN, results_Track_VAN, results_Track_DAN_VAN] = getResults(Condition,Net1,Net2);

Net1 = 'VAN';
Net2 = 'VIS';
[results_Track_VAN, results_Track_VIS, results_Track_VAN_VIS] = getResults(Condition,Net1,Net2);

Net1 = 'DAN';
Net2 = 'VIS';
[results_Track_DAN, results_Track_VIS, results_Track_DAN_VIS] = getResults(Condition,Net1,Net2);

plotCondition(Condition,results_Track_DAN,results_Track_VAN,results_Track_VIS,results_Track_DAN_VAN,results_Track_VAN_VIS,results_Track_DAN_VIS);

end

function [results_Net1, results_Net2, results_Net1_Net2] = getResults(Condition,Net1,Net2)

load(strcat('LHR','-','SUBJ2','-','Granger-Functional-Seeds','-',Condition,'-',Net1,'-',Net2,'.mat'));

nVoxels_net1 = length(all_voxels_inet)
nVoxels_net2 = length(all_voxels_iinet)

size(Nets_Gpval)

pcriterion = 0.01;

Gpval_net1 = Nets_Gpval(:,1:nVoxels_net1,1:nVoxels_net1);
Gpval_net2 = Nets_Gpval(:,(nVoxels_net1+1):(nVoxels_net1+nVoxels_net2),(nVoxels_net1+1):(nVoxels_net1+nVoxels_net2));
Gpval_net1_net2 = Nets_Gpval(:,(nVoxels_net1+1):(nVoxels_net1+nVoxels_net2),1:nVoxels_net1);
Gpval_net2_net1 = Nets_Gpval(:,1:nVoxels_net1,(nVoxels_net1+1):(nVoxels_net1+nVoxels_net2));

results_Net1 = length(find(Gpval_net1 < pcriterion)) / length(Gpval_net1(:));
results_Net2 = length(find(Gpval_net2 < pcriterion)) / length(Gpval_net2(:));

results_Net1_Net2 = ( length(find(Gpval_net1_net2 < pcriterion)) + length(find(Gpval_net2_net1 < pcriterion)) ) / ( length(Gpval_net1_net2(:)) + length(Gpval_net2_net1(:)) );

end

function plotCondition(Condition,DAN,VAN,VIS,DAN_VAN,VAN_VIS,DAN_VIS)

matrix(1,1) = DAN;
matrix(2,1) = DAN_VAN;
matrix(3,1) = DAN_VIS;

matrix(1,2) = DAN_VAN;
matrix(2,2) = VAN;
matrix(3,2) = VAN_VIS;

matrix(1,3) = DAN_VIS;
matrix(2,3) = VAN_VIS;
matrix(3,3) = VIS;

figure;

imagesc(matrix);

colorbar;

hold on

end
