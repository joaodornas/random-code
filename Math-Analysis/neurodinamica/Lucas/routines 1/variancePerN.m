variance_Pref_SF=zeros(size(Pref_SF,1)-1,1);
for i=1:size(Pref_SF,1)-1
    variance_Pref_SF(i)=var(Pref_SF(1:i+1));
end

figure;
subplot(2,3,1)
plot(variance_Pref_SF);
title('variance Pref-SF');

variance_Pref_TF=zeros(size(Pref_TF,1)-1,1);
for i=1:size(Pref_TF,1)-1
    variance_Pref_TF(i)=var(Pref_TF(1:i+1));
end

subplot(2,3,2)
plot(variance_Pref_TF);
title('variance Pref-TF');

variance_SF_BW=zeros(size(SF_BW,1)-1,1);
for i=1:size(SF_BW,1)-1
    variance_SF_BW(i)=var(SF_BW(1:i+1));
end

subplot(2,3,3)
plot(variance_SF_BW);
title('variance SF-BW');

variance_SF_resolution=zeros(size(SF_resolution,1)-1,1);
for i=1:size(SF_resolution,1)-1
    variance_SF_resolution(i)=var(SF_resolution(1:i+1));
end

subplot(2,3,4)
plot(variance_SF_resolution);
title('variance SF-resolution');

variance_TF_BW=zeros(size(TF_BW,1)-1,1);
for i=1:size(TF_BW,1)-1
    variance_TF_BW(i)=var(TF_BW(1:i+1));
end

subplot(2,3,5)
plot(variance_TF_BW);
title('variance TF-BW');

variance_TF_resolution=zeros(size(TF_resolution,1)-1,1);
for i=1:size(TF_resolution,1)-1
    variance_TF_resolution(i)=var(TF_resolution(1:i+1));
end

subplot(2,3,6)
plot(variance_TF_resolution);
title('variance TF-resolution');