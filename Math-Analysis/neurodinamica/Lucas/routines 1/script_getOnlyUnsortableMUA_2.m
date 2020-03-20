BaselineSEM_a=BaselineSEM;
STTCMeanResponse_a=STTCMeanResponse;
baseline_a=baseline;
twoSEM_a=twoSEM;
Sf_cpd=2.^Sf_octave;

load('multiple_stt021a02_1b','STTCMeanResponse');
STTCMeanResponse_b=STTCMeanResponse;
load('SF_stt021a02_1b','baseline');
baseline_b=baseline;
load('multiple_stt021a02_1c','STTCMeanResponse');
STTCMeanResponse_c=STTCMeanResponse;
load('SF_stt021a02_1c','baseline');
baseline_c=baseline;
load('multiple_stt021a02_1d','STTCMeanResponse');
STTCMeanResponse_d=STTCMeanResponse;
load('SF_stt021a02_1d','baseline');
baseline_d=baseline;

[ibest_sf,ibest_tf]=find(STTCMeanResponse_a==max(max(STTCMeanResponse_a)));
SF_curve=STTCMeanResponse_a(:,ibest_tf)';
TF_curve=STTCMeanResponse_a(ibest_sf,:);

ncond=size(Sf_octave,2)*size(Tf_octave,2);

vector_mua=zeros(ncond,1);
vector_sua1=zeros(ncond,1);
vector_sua2=zeros(ncond,1);
vector_sua3=zeros(ncond,1);

for i=1:ncond
    vector_mua(i)=STTCMeanResponse_a(i)/max(max(STTCMeanResponse_a));
    vector_sua1(i)=STTCMeanResponse_b(i)/max(max(STTCMeanResponse_b));
    vector_sua2(i)=STTCMeanResponse_c(i)/max(max(STTCMeanResponse_c));
    vector_sua3(i)=STTCMeanResponse_d(i)/max(max(STTCMeanResponse_d));
end

[r1,p1]=corr(vector_mua,vector_sua1,'type','Spearman','rows','pairwise')
[r2,p2]=corr(vector_mua,vector_sua2,'type','Spearman','rows','pairwise')
[r3,p3]=corr(vector_mua,vector_sua3,'type','Spearman','rows','pairwise')

figure;
contourf(STTCMeanResponse_a');
figure;
contourf(STTCMeanResponse_d');