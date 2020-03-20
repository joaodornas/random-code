BaselineSEM_mina=BaselineSEM;
STTCMeanResponse_mina=STTCMeanResponse;
baseline_mina=baseline;
twoSEM_mina=twoSEM;
Sf_cpd=2.^Sf_octave;

load('multiple_stt006b02_1a','STTCMeanResponse');
STTCMeanResponse_a=STTCMeanResponse;
load('SF_stt006b02_1a','baseline');
baseline_a=baseline;

baseline_min_AandB=baseline_mina-baseline_b
STTCMeanResponse_min_AandB=((STTCMeanResponse_mina+baseline_mina)-(STTCMeanResponse_b+baseline_b))-baseline_min_AandB;

[ibest_sf,ibest_tf]=find(STTCMeanResponse_min_AandB==max(max(STTCMeanResponse_min_AandB)));
SF_curve=STTCMeanResponse_min_AandB(:,ibest_tf)';
TF_curve=STTCMeanResponse_min_AandB(ibest_sf,:);

ncond=size(Sf_octave,2)*size(Tf_octave,2);

vector_mua=zeros(ncond,1);
vector_sua=zeros(ncond,1);
vector_sua2=zeros(ncond,1);

for i=1:ncond
    %vector_mua(i)=STTCMeanResponse_min_AandB(i)/max(max(STTCMeanResponse_min_AandB));
    vector_sua2(i)=STTCMeanResponse_a(i)/max(max(STTCMeanResponse_a));
end

[r2,p2]=corr(vector_mua,vector_sua2,'type','Spearman','rows','pairwise')

figure;
contourf(STTCMeanResponse_min_AandB');
figure;
contourf(STTCMeanResponse_b');