clear all
load '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Others\population_data'

mean_Pref_SF=zeros(size(Pref_SF,1)-1,1);
for i=1:size(Pref_SF,1)-1
    mean_Pref_SF(i)=mean(Pref_SF(1:i+1));
end

figure;
subplot(2,3,1)
plot(mean_Pref_SF);
title('mean Pref-SF');

mean_Pref_TF=zeros(size(Pref_TF,1)-1,1);
for i=1:size(Pref_TF,1)-1
    mean_Pref_TF(i)=mean(Pref_TF(1:i+1));
end

subplot(2,3,2)
plot(mean_Pref_TF);
title('mean Pref-TF');

mean_SF_BW=zeros(size(SF_BW,1)-1,1);
for i=1:size(SF_BW,1)-1
    mean_SF_BW(i)=mean(SF_BW(1:i+1));
end

subplot(2,3,3)
plot(mean_SF_BW);
title('mean SF-BW');

mean_SF_resolution=zeros(size(SF_resolution,1)-1,1);
for i=1:size(SF_resolution,1)-1
    mean_SF_resolution(i)=mean(SF_resolution(1:i+1));
end

subplot(2,3,4)
plot(mean_SF_resolution);
title('mean SF-resolution');

mean_TF_BW=zeros(size(TF_BW,1)-1,1);
for i=1:size(TF_BW,1)-1
    mean_TF_BW(i)=mean(TF_BW(1:i+1));
end

subplot(2,3,5)
plot(mean_TF_BW);
title('mean TF-BW');

mean_TF_resolution=zeros(size(TF_resolution,1)-1,1);
for i=1:size(TF_resolution,1)-1
    mean_TF_resolution(i)=mean(TF_resolution(1:i+1));
end

subplot(2,3,6)
plot(mean_TF_resolution);
title('mean TF-resolution');

save '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Others\meanPerN'