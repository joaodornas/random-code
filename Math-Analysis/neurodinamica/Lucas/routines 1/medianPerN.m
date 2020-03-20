clear all
load '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Others\population_data'

median_Pref_SF=zeros(size(Pref_SF,1)-1,1);
for i=1:size(Pref_SF,1)-1
    median_Pref_SF(i)=median(Pref_SF(1:i+1));
end

figure;
subplot(2,3,1)
plot(median_Pref_SF);
title('median Pref-SF');

median_Pref_TF=zeros(size(Pref_TF,1)-1,1);
for i=1:size(Pref_TF,1)-1
    median_Pref_TF(i)=median(Pref_TF(1:i+1));
end

subplot(2,3,2)
plot(median_Pref_TF);
title('median Pref-TF');

median_SF_BW=zeros(size(SF_BW,1)-1,1);
for i=1:size(SF_BW,1)-1
    median_SF_BW(i)=median(SF_BW(1:i+1));
end

subplot(2,3,3)
plot(median_SF_BW);
title('median SF-BW');

median_SF_resolution=zeros(size(SF_resolution,1)-1,1);
for i=1:size(SF_resolution,1)-1
    median_SF_resolution(i)=median(SF_resolution(1:i+1));
end

subplot(2,3,4)
plot(median_SF_resolution);
title('median SF-resolution');

median_TF_BW=zeros(size(TF_BW,1)-1,1);
for i=1:size(TF_BW,1)-1
    median_TF_BW(i)=median(TF_BW(1:i+1));
end

subplot(2,3,5)
plot(median_TF_BW);
title('median TF-BW');

median_TF_resolution=zeros(size(TF_resolution,1)-1,1);
for i=1:size(TF_resolution,1)-1
    median_TF_resolution(i)=median(TF_resolution(1:i+1));
end

subplot(2,3,6)
plot(median_TF_resolution);
title('median TF-resolution');

save '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Others\medianPerN'