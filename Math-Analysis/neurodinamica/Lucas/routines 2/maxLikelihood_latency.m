function [lat,cutoff]=maxLikelihood_latency(psth,parameters)

%[lat,cutoff]=maxLikelihood_latency(psth,parameters)
%calculates latency (lat) with the maximum likelihood estimation method
%described in Friedman & Priebe (1998) J Neurosci Meth 83(2):185-194
%
%inputs: psth: vector containing psth values (typically bin = 1 ms)
%parameters: standard matlab structure used in dinamo, containing the
%following fields: 
%baseline_duration: duration of baseline in ms, default 1000
%analysis_period: period of the psth to be analyzed in ms, default 300
%analysis_startTime: time after stimulus onset to start latency calculation
%(in ms), default 0
%
%outputs: lat is latency in ms, cutoff is the estimated point where
%response becomes stationary (in ms), and is used in the maximum likelihood
%estimation. See reference for details
%
%by Lucas Pinto in April 2009


if nargin==1
    baseline_duration=1000;
    analysis_period=200;
    analysis_startTime=20;
elseif nargin==2
    baseline_duration=parameters.baseline_duration;
    analysis_period=parameters.analysis_period;
    analysis_startTime=parameters.analysis_startTime;
elseif nargin < 1
    error 'i need inputs!'
elseif nargin > 2
    error 'too many inputs - see help'
end

%truncate psth
trunc_psth=psth(baseline_duration+1:baseline_duration+analysis_period);

%find cutoff that maximizes the difference between the slopes of two line
%segments corresponding to baseline and response periods in the cumulative psth
cum_psth=cumsum(trunc_psth);

%initialize vectors
slopediff1=[];
% slopediff2=[];
candidate_lat=[];
candidate_cutoff=[];

for k=analysis_startTime+4:size(cum_psth,2)-2
    for t=analysis_startTime+2:k-2
        %slope of the line that goes from start time to latency  
        line1=cum_psth(1:t);
        p1=polyfit(1:size(line1,2),line1,1);
        slope1=p1(1);
        
        %slope of the line that goes from latency to cutoff 
        line2=cum_psth(t+1:k);
        p2=polyfit(1:size(line2,2),line2,1);
        slope2=p2(1);
        
%         %slope of the line that goes from cutoff to end of analysis period
%         line3=cum_psth(k+1:size(cum_psth,2));
%         p3=polyfit(1:size(line3,2),line3,1);
%         slope3=p3(1);
        
        slopediff1=[slopediff1 slope2-slope1];
%         slopediff2=[slopediff2 slope2-slope3];
        candidate_lat=[candidate_lat t];
        candidate_cutoff=[candidate_cutoff k];
    end
end

maxslopediff1_i=find(slopediff1==max(slopediff1),1,'first');
% maxslopediff2_i=find(slopediff2==max(slopediff2));
% lat_fromSlope=candidate_lat(maxslopediff1_i(1));
% cutoff=candidate_cutoff(maxslopediff2_i(1));
cutoff=candidate_cutoff(maxslopediff1_i);

%manipulate psth to avoid error in factorial calculation
maxvalue=max(trunc_psth);
if maxvalue >=100
    trunc_psth=fix(trunc_psth./100);
elseif maxvalue<=1
    trunc_psth=fix(trunc_psth.*10);
elseif maxvalue>1 && maxvalue<=10
    trunc_psth=fix(trunc_psth);
elseif maxvalue>10 && maxvalue<=100
    trunc_psth=fix(trunc_psth./10);
end

%find latency value (l) that maximizes log-likelihood function written below
likelihood=zeros(1,cutoff-1-analysis_startTime);

for l=analysis_startTime:cutoff-1

    baseline_rate = sum(trunc_psth(1:l))/(l+1);
    response_rate = sum(trunc_psth(l+1:cutoff))/(cutoff-l);

    likelihood(l-analysis_startTime+1) = - baseline_rate*(l+1) + log(baseline_rate)*sum(trunc_psth(1:l)) - sum(log(factorial(trunc_psth(1:l)))) - response_rate*(cutoff-l) + log(response_rate).*sum(trunc_psth(l+1:cutoff)) - sum(log(factorial(trunc_psth(l+1:cutoff))));

end

lats=find(likelihood==max(likelihood),1,'first');
lat=lats+analysis_startTime-1;