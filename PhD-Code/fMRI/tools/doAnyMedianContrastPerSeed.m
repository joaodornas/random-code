function [contrast_attention_z, contrast_stimulus_z, p_attention, p_stimulus, h_attention, h_stimulus] = doAnyMedianContrastPerSeed(my_samples,level,label)

nTotalClusters = 758;
pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

for iIg=1:N
    
   disp(strcat('Ig:',int2str(iIg)));
      
   p_attention(iIg) = ranksum(squeeze(my_samples(:,idx_Attention,iIg)),squeeze(my_samples(:,idx_PassiveViewing,iIg)));
   
   if median(squeeze(my_samples(:,idx_Attention,iIg))) > median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); h_attention(iIg) = -1; end
   if median(squeeze(my_samples(:,idx_Attention,iIg))) < median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); h_attention(iIg) = 1; end
  
   p_stimulus(iIg) = ranksum(squeeze(my_samples(:,idx_PassiveViewing,iIg)),squeeze(my_samples(:,idx_RestingState,iIg)));
   
   if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) > median(squeeze(my_samples(:,idx_RestingState,iIg))); h_stimulus(iIg) = -1; end
   if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) < median(squeeze(my_samples(:,idx_RestingState,iIg))); h_stimulus(iIg) = 1; end

end

contrast_attention = zeros(1,N);
contrast_stimulus = zeros(1,N);

contrast_attention_z = zeros(1,N);
contrast_stimulus_z = zeros(1,N);

for iIg=1:N
      
    contrast_attention(iIg) = p_attention(iIg);
    contrast_stimulus(iIg) = p_stimulus(iIg);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iIg=1:N
     
    if contrast_attention(iIg) ~= 0; contrast_attention_z(iIg) = norminv(contrast_attention(iIg)); end
    if contrast_stimulus(iIg) ~= 0; contrast_stimulus_z(iIg) = norminv(contrast_stimulus(iIg)); end

    contrast_attention_z(iIg) = h_attention(iIg) * contrast_attention_z(iIg);
    contrast_stimulus_z(iIg) = h_stimulus(iIg) * contrast_stimulus_z(iIg);
    
end

% save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'contrast_stimulus_z','contrast_attention_z');

end