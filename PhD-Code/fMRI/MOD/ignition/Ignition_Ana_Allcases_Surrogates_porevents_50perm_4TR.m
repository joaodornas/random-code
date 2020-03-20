clear all;
dir = 'E:\BCN_ANALYSIS\Drugs\PSY_LSD_SLEEP\LSD';  %put your directoy
cd(dir)
%load SC90_LRLR.mat;
% C=SC_new;
% %%%%
% C=C/max(max(C))*0.2;
% Cbin=C>0;

load LSD_aal.mat;   %cell with nsubjectsXconditions


for cond= 1:6
 perm = 50;
 
 for permut = 1:perm

msg = sprintf('CASE %d permutation  %d', cond, permut) 
 

%adapt parameters TO YOUR DATA
TR= 2; %sampling interval(TR)  
Tmax=217; %timepoints
N=90;%758;%90;  regions
NSUB=15; %subjects
Isubdiag = find(tril(ones(N),-1));



 CASE=cond;




TC1=tc_aal{1,CASE};
TC2=tc_aal{2,CASE};
TC3=tc_aal{3,CASE};
TC4=tc_aal{4,CASE};
TC5=tc_aal{5,CASE};
TC6=tc_aal{6,CASE};
TC7=tc_aal{7,CASE};
TC8=tc_aal{8,CASE};
TC9=tc_aal{9,CASE};
TC10=tc_aal{10,CASE};
TC11=tc_aal{11,CASE};
TC12=tc_aal{12,CASE};
TC13=tc_aal{13,CASE};
TC14=tc_aal{14,CASE};
TC15=tc_aal{15,CASE};  %as many subjects as in your data


%%%%%%%%%%%%%
nevents=zeros(1,N);


 FC=zeros(NSUB,N,N);

 flp = .04;           % lowpass frequency of filter
 fhi = .07;           % highpass
 delt = TR;            % sampling interval
 k=2;                  % 2nd order butterworth filter
 fnq=1/(2*delt);       % Nyquist frequency
 Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
 [bfilt2,afilt2]=butter(k,Wn);   % construct the filter
   
 ttotal=1;

 
 
 
     
  for nsub=1:NSUB
    xs=eval(sprintf('TC%d', nsub));
    Tmax=size(xs,2); 
    T=10:Tmax-10;  
   for seed=1:N
    x=demean(detrend(xs(seed,:)));
    timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
    Xanalytic = hilbert(demean(timeseriedata(seed,:)));
    tise = detrend(demean(timeseriedata(seed,T)));
    ev1=tise>std(tise)+mean(tise);
    ev2=[0 ev1(1:end-1)];
    events(seed,:)=(ev1-ev2)>0; 
    events_perm(seed, :) = events(seed,randperm(size(events(seed,:),2)));
    Phases(seed,:) = angle(Xanalytic);
  end

 %%% integration
 
%  T=1:1:size(Phases,2);

 for t=T
  for i=1:N
    for j=1:N
   %  phasematrix(i,j)=exp(-3*adif(Phases(i,t),Phases(j,t)));
   %  phasematrix(i,j)=events(i,t-9)*events(j,t-9);
     phasematrix_perm(i,j)=events_perm(i,t-9)*events_perm(j,t-9);
    end
  end 
 % cc=phasematrix;%*Cbin;
  cc=phasematrix_perm;%*Cbin;
  cc=cc-eye(N);
  
    [comps csize]=get_components(cc);
    integ(t-9)=max(csize)/N;

%   pp=1;
%   PR=0:0.01:0.99;
%   for p=PR
%    A=abs(cc)>p;
%    [comps csize]=get_components(A);
%    cs(pp)=max(csize);
%    pp=pp+1;
%   end
%   integ(t-9)=sum(cs)*0.01/N;
  
%   integt(ttotal)=sum(cs)*0.01/N;
%   [M modularity(t-9)]=community_louvain(phasematrix);
%  ttotal=ttotal+1;
 end
 
 %%%
   
 %%%% event trigger
 nevents2=zeros(1,N);

  for seed=1:N
   flag=0;
   for t=T
    %if events(seed,t-9)==1 %%&& flag==0
    if events_perm(seed,t-9)==1 && flag==0
     flag=1;
     nevents(seed)=nevents(seed)+1;
     nevents2(seed)=nevents2(seed)+1;
    end
    if flag>0
     IntegStim(seed,flag,nevents(seed))=integ(t-9);
     IntegStim2(seed,flag,nevents2(seed))=integ(t-9);
%      IntegStimQ(seed,flag,nevents(seed))=modularity(t-9);
     flag=flag+1;
    end
    if flag==5
      flag=0;
     end
   end
  end
 
%   for seed=1:N
%    mevokedinteg2(seed)=max(mean(squeeze(IntegStim2(seed,:,1:nevents2(seed))),2));
%   end
%   mignitionon2(nsub)=mean(mevokedinteg2);
%   stdignitionon2(nsub)=std(mevokedinteg2);

end

% mignitionon=mean(mignitionon2);
% stdignitionon=mean(stdignitionon2);

%minteg=mean(integt);
for seed=1:N
 %mevokedinteg(seed,:, permut)=mean(squeeze(IntegStim(seed,:,1:nevents(seed))),2);
 mev_byperm(seed, permut)=mean(squeeze(max(squeeze(IntegStim(seed,:,1:nevents(seed))),[],1)));
 stdev_byperm(seed, permut)=std(squeeze(max(squeeze(IntegStim(seed,:,1:nevents(seed))),[],1)));
%  mevokedQ(seed,:)=mean(squeeze(IntegStimQ(seed,:,1:nevents(seed))),2);
end

clearvars -except tc_aal Cbin cond permut  mev_byperm stdev_byperm
 end

  
%for seed=1:N
 %mevokedinteg(seed,:)= mean(mevokedinteg,3);
 mev=mean(mev_byperm, 2);
 stdev= mean (stdev_byperm,2);
%  mevokedQ(seed,:)=mean(squeeze(IntegStimQ(seed,:,1:nevents(seed))),2);
%end


 
save (sprintf('4TR_%dpermut_Ignition_surroRandEv_CASE%d.mat', permut, cond), 'mev', 'stdev');
clearvars -except tc_aal %Cbin
end  
  