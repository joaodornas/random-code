function [contrast_attention_z, contrast_stimulus_z, p_attention, p_stimulus, h_attention, h_stimulus] = doAnyMedianContrastPerParam(my_samples,level,label)

disp(strcat('isINF:',num2str(length(find(isinf(abs(my_samples)))) / length(my_samples(:)))));
my_samples(find(isinf(abs(my_samples)))) = 0;

nTotalClusters = 758;
pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

h_attention = 0;
h_stimulus = 0;

p_attention = ranksum(squeeze(my_samples(:,idx_Attention)),squeeze(my_samples(:,idx_PassiveViewing)));

if median(squeeze(my_samples(:,idx_Attention))) > median(squeeze(my_samples(:,idx_PassiveViewing))); h_attention = -1; end
if median(squeeze(my_samples(:,idx_Attention))) < median(squeeze(my_samples(:,idx_PassiveViewing))); h_attention = 1; end

p_stimulus = ranksum(squeeze(my_samples(:,idx_PassiveViewing)),squeeze(my_samples(:,idx_RestingState)));

if median(squeeze(my_samples(:,idx_PassiveViewing))) > median(squeeze(my_samples(:,idx_RestingState))); h_stimulus = -1; end
if median(squeeze(my_samples(:,idx_PassiveViewing))) < median(squeeze(my_samples(:,idx_RestingState))); h_stimulus = 1; end

contrast_attention = p_attention;
contrast_stimulus = p_stimulus;

if contrast_attention ~= 0; contrast_attention_z = norminv(contrast_attention); end
if contrast_stimulus ~= 0; contrast_stimulus_z = norminv(contrast_stimulus); end

contrast_attention_z = h_attention * contrast_attention_z;
contrast_stimulus_z = h_stimulus * contrast_stimulus_z;

% save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'contrast_stimulus_z','contrast_attention_z');

end
