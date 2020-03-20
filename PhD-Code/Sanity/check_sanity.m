function check_sanity

clear all;
close all


%% run granger causality on portions of sanity-check data


% retrieve values Ci and Di to file 

load( 'sanity.mat', 'Ci', 'Di', 'Tevent' );

% time axis

Tstart = 1;
Tend   = length(Ci);
Ti     = 1:Tend;

% event number

Nevent = length(Tevent);

% number,size, density of samples

Nsample = 50;

Nsize   = 2000;

Density = 5;  % one in five 


% loop over samples

maxlag      = 5;


GC_CtoD = nan( Nsample, 2 );
GC_DtoC = nan( Nsample, 2 );

for n=1:Nsample
    
    % choose sample
    
    Sstart = round( rand * ( Tend - Nsize * Density ) );
    
    Send   = Sstart + (Nsize-1) * Density;
    
    idxS   = Sstart:Density:Send;
    [F,Cv] = granger_cause(Di(idxS),Ci(idxS),0.05,maxlag);
    GC_CtoD(n,:) =  [F,Cv];  
    [F,Cv] = granger_cause(Ci(idxS),Di(idxS),0.05,maxlag);   
    GC_DtoC(n,:) =  [F,Cv];
    
end

mxF  = max([GC_CtoD(:,1)' GC_DtoC(:,1)']);

mxCV = max([GC_CtoD(:,2)' GC_DtoC(:,2)']);

mnF  = min([GC_CtoD(:,1)' GC_DtoC(:,1)']);

mnCV = min([GC_CtoD(:,2)' GC_DtoC(:,2)']);

avF  = mean([GC_CtoD(:,1)' GC_DtoC(:,1)']);

avCV = mean([GC_CtoD(:,2)' GC_DtoC(:,2)']);

Frange = [mnF avF mxF]

CVrange = [mnCV avCV mxCV]

% plot results


figure;

X = linspace(1,10,41);

[Ncd,X] = hist( log(GC_CtoD(:,1)), X );

[Ndc,X] = hist( log(GC_DtoC(:,1)), X );

hold on;



h = bar( X, Ndc );

set(h, 'FaceColor', [0 0 1]);

h = bar( X, Ncd );

set(h, 'FaceColor', [1 0 0]);

plot( log(avCV)*[1 1], [0 40], 'k--', 'LineWidth', 2 );

hold off;


% % loop over transitions
% 
% noiseflag = 1;
% alpha     = 0.5;  % noise sigma
% 
% maxlag    = 5;
% 
% for it = 1:Ntrans
%     
%    itr = itrans(it); 
%     
%    % loop over granger windows
%    
%    for ik = 1:Noffset
%        
%        kos = kwoffset(ik);
%        
%        % crop windows
%        
%        kk = itr+kos+iwsample;  % index to window entries
%        
%        kE1 = dE1(kk);
%        kE2 = dE2(kk);
%        kD1 = dD1(kk);
%        kD2 = dD2(kk);
%        kE1m2 = dE1m2(kk);
%        
%        % subtract mean
%        
%        kE1   = kE1 - mean(kE1);
%        kE2   = kE2 - mean(kE2);
%        kD1   = kD1 - mean(kD1);
%        kD2   = kD2 - mean(kD2);
%        kE1m2 = kE1m2 - mean(kE1m2);
%        
%        % divide by standard deviation (if nonzero)
%        
%        if std(kD1) > 0   % effect vector has variance
%            kD1   = kD1 / std(kD1);
%            if noiseflag
%                kD1 = kD1 + normrnd( 0, alpha, size(kk) );
%            end
%            
%            
% %            if std(kE1) > 0  % cause vector has variance
% %                kE1   = kE1 / std(kE1);
% %                if noiseflag
% %                    kE1 = kE1 + normrnd( 0, alpha, size(kk) );
% %                end
% %                E1gcD1(it,ik,:)   = granger_cause(kD1,kE1,0.05,maxlag);
% %            end
% %            
% %            if std(kE2) > 0
% %                kE2   = kE2 / std(kE2);
% %                if noiseflag
% %                    kE1 = kE1 + normrnd( 0, alpha, size(kk) );
% %                end
% %                E2gcD1(it,ik,:)   = granger_cause(kD1,kE2,0.05,maxlag);
% %            end
%            
%            if std(kD2) > 0
%                kD2   = kD2 / std(kD2);
%                if noiseflag
%                    kD2 = kD2 + normrnd( 0, alpha, size(kk) );
%                end
%                D2gcD1(it,ik,:)   = granger_cause(kD1,kD2,0.05,maxlag);
%            end
%            
%            if std(kE1m2) > 0
%                kE1m2   = kE1m2 / std(kE1m2);
%                if noiseflag
%                    kE1m2 = kE1m2 + normrnd( 0, alpha, size(kk) );
%                end
%                E1m2gcD1(it,ik,:) = granger_cause(kD1,kE1m2,0.05,maxlag);
%            end
%            
%        end
%        
%    end
%    
%    if mod(it,100)==0
%        it
%    end
%    
% end
% 
% %% average results over transitions
% 
% % average F-value
% 
% aveF_E1      = nanmean( E1gcD1( :, :, 1), 1 );
% aveF_E2      = nanmean( E2gcD1( :, :, 1), 1 );
% aveF_D2      = nanmean( D2gcD1( :, :, 1), 1 );
% aveF_E1m2    = nanmean( E1m2gcD1( :, :, 1), 1 );
% 
% % average rejected null hypothesis (confirmed GC)
% 
% aveGC_E1     = nanmean( E1gcD1(:,:,1) > E1gcD1(:,:,2), 1 );
% aveGC_E2     = nanmean( E2gcD1(:,:,1) > E2gcD1(:,:,2), 1 );
% aveGC_D2     = nanmean( D2gcD1(:,:,1) > D2gcD1(:,:,2), 1 );
% aveGC_E1m2   = nanmean( E1m2gcD1(:,:,1) > E1m2gcD1(:,:,2), 1 );
% 
% 
% 
% 
% 
% return;