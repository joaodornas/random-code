function h_v2


%% Load Data Sets

dataLength = 900000;

[E, D, Trans, dataLabel] = getData('Sanity',dataLength);

%[E, D, Trans, dataLabel] = getData('Robin',dataLength);

%[E, D, Trans, dataLabel] = getData('Random',dataLength);

%[E, D, Trans, dataLabel] = getData('Senoide',dataLength);

%% Get Segments Around Transitions with Hilbert

%[E_seg, D_seg] = getAndPlotSegmentsAroundTransitions(E,D,Trans,dataLabel);


%% Get Hilbert Transform - Full Time Series

[hEa, hEw, hEwt, hEy] = hilbertTransform(E);
[hDa, hDw, hDwt, hDy] = hilbertTransform(D);

%% Pre-Process the Data

%[hEw, hDw] = absDegreePhases(hEw,hDw);

%[hEw, hDw] = getAndPlotWalshHadamard(hEw,hDw,Trans,2000);

%[hEw, hDw] = getDetrend(hEw,hDw,Trans);

%[hEw, hDw] = smoothAll(hEw,hDw,100); 

%[hEw, hDw] = unwrapAll(hEw,hDw);

%[hEw, hDw] = unwrapAllAlgorithm(hEw,hDw);

%[hEw, hDw] = setLocalMeanToZero(hEw,hDw);

%[hEw, hDw, Trans] = onlyASegment(hEw,hDw,Trans,2000);


%% Plot Hilbert Transform Information

%plotHilbert(E,D,hDa,hEa,hDw,hEw,hDwt,hEwt,Trans,dataLabel);

%% Plot Phase Distribution At Transitions

%plotPhaseDistributionAtTransitions(hEw,hDw,Trans,dataLabel);

%% Get and Plot Auto and Cross Correlations

%getAndPlotAutoAndCrossCorrelations(E,D,hEw,hDw,dataLabel);

%% Get and Plot Conditional Expectations

getAndPlotConditionalExpectations(hEw,hDw,Trans,dataLabel);

%% Granger Causality

%getGranger(hEw,hDw,dataLabel);

%% Get and Plot Relative Phases

%getAndPlotRelativePhases(hEw,hDw,Trans,dataLabel);

%% Get and Plot Orbits

%getAndPlotOrbitsAndRelativePhases(hEw,hDw,Trans,dataLabel);

%% Get and Plot Relative Phase Distributions At Transitions

%plotRelativePhaseDistributionAtTrantisionsWithTau(hEw,hDw,Trans,dataLabel);

%% Get and Plot Detrend

%getAndPlotDetrend(hEw,hDw,Trans,dataLabel)

%% Plot Phases in Time

%plotPhasesInTime(hEw,hDw,Trans,dataLabel);

%% Plot Phase Distribution at Tau

%plotPhaseDistributionAtTau(hEw,hDw,Trans,dataLabel);

%% Plot Orbits Phases Relative to Transition

%Trans = Trans(1);

%plotOrbitsPhasesRelativeToTransitions(hEw,hDw,Trans,dataLabel);

%% Plot Phases Differences Distribution at Transitions

%plotPhaseDifferenceDistributionAtTransition(hEw,hDw,Trans,dataLabel)

%% Plot Fourier Frequencies Amplitudes

%plotFourierFrequenciesAmplitudes(hEw,hDw,dataLabel);

%% Get And Plot Condition Phase Distribution

%getAndPlotConditionalPhaseDistribution(hEw,hDw,Trans,dataLabel);
%getAndPlotConditionalPhaseDistribution(hDw,hEw,Trans,dataLabel);

%% Get and Plot Correlations of Phases Shifted in Time by Tau

%[tau_e, tau_d, max_tau] = plotCorrelationsForManyTau(hEw,hDw,dataLabel);

%% Get Transitions Segments

%[hEwseg, hDwseg] = getTransitionsSegments(hEw,hDw,Trans);


%% Get and Plot Causation At a Tau Offset in Segments

%getAndPlotCausationAtTauOffsetInSegments(hEwseg,hDwseg,400);


function [hEw, hDw, Trans] = onlyASegment(hEw,hDw,Trans,size)
        
    nTrans = length(Trans);
    
    rand_trans = Trans(randperm(nTrans)); 
    
    idx_trans = rand_trans(1);
    
    offset = -(size/2):(size/2);
    
    hEw = hEw(offset + idx_trans);
    hDw = hDw(offset + idx_trans);
    
    Trans = [];
    Trans(1) = size/2 + 1;
        
end

function getAndPlotCausationAtTauOffsetInSegments(hEwseg,hDwseg,tau)
        
        o = 1;
    
        for i=(tau+1):length(hEwseg(1).seg)
            
           for j=1:length(hEwseg)
               
                phases_e(j) = hEwseg(j).seg(i - tau);
                phases_d(j) = hDwseg(j).seg(i);
  
           end
           
           [unique_e, nUnique_e] = count_unique(phases_e);
           
           for k=1:length(unique_e)
               
               idx_unique_e = find(phases_e == unique_e(k));
               
               selected_phases_d = phases_d(idx_unique_e);
               
               [unique_d, nUnique_d] = count_unique(selected_phases_d);
               
               prob_density = nUnique_d./sum(nUnique_d);
               
               [Chance(k), H, n] = getEntropy(prob_density);
               
               Chance(k) = Chance(k)*nUnique_e(k);
          
           end
 
           Causality(o) = sum(Chance)/sum(unique_e);
           
           o = o + 1;
           
        end
        
        idx_trans = floor(length(hEwseg(1).seg)/2);
        
        figure
        plot(1:length(Causality),Causality,'b');
        title(strcat(dataLabel,':Causation'));
        hold on
        plot(idx_trans,max(Causality),'ro');
        hold on
        plot([idx_trans idx_trans],[0 max(Causality)],'r');
    
end

function [hEwseg, hDwseg] = getTransitionsSegments(hEw,hDw,Trans)
        
    half_size = 1000;
    
    offset = -half_size:half_size;
    
    j = 0;
    for i=1:length(Trans)
        
        if ((Trans(i) - half_size) > 0) && ((Trans(i) + half_size) < length(hEw))
            
            j = j + 1;
            
            hEwseg(j).seg = hEw(offset + Trans(i));
            hDwseg(j).seg = hDw(offset + Trans(i));
        
        end
  
    end
     
end

function [Causality, H, N] = getCausality(Chance)
       
        H = 0;
 
        if length(Chance) == 1
            
            H = 0;
            Causality = Chance(1);
            N = 1;
            
        else
            
            for i=1:length(Chance)

                H = -1 * Chance(i) * log2( Chance(i) ) + H;

            end

           N = length(Chance);

           Causality = 1 - H/log2(N);

        end
        
end

function [Chance, H, N] = getEntropy(p)
       
        H = 0;
 
        if length(p) == 1
            
            H = 0;
            Chance = p(1);
            N = 1;
            
        else
            
            for i=1:length(p)

                H = -1 * p(i) * log2( p(i) ) + H;

            end

           N = length(p);

           Chance = 1 - H/log2(N);

        end
        
end

function [tau_e, tau_d, max_tau] = plotCorrelationsForManyTau(hEw,hDw,dataLabel)
    
    max_tau = 2000;
    alpha = 0.05;
    
    for tau=1:max_tau
    
        phases_e = hEw((max_tau + 1):end);
        phases_d_s = hDw((max_tau - tau + 1):(end - tau));
        
        phases_d = hDw((max_tau + 1):end);
        phases_e_s = hEw((max_tau - tau+ 1):(end - tau));
        
        phases_e = setToColumnVector(phases_e);
        phases_d_s = setToColumnVector(phases_d_s);
        
        phases_d = setToColumnVector(phases_d);
        phases_e_s = setToColumnVector(phases_e_s);
        
        [r_d_s(tau), p_d_s(tau)] = corr(phases_e,phases_d_s);
        [r_e_s(tau), p_e_s(tau)] = corr(phases_d,phases_e_s);
        
    end
    
    idx_max_r_d_s = find(r_d_s == max(r_d_s));
    idx_max_r_e_s = find(r_e_s == max(r_e_s));
    
    idx_d_s_sig = p_d_s < alpha;
    idx_e_s_sig = p_e_s < alpha;
    
    figure
    plot(1:max_tau,r_d_s,'b');
    hold on
    plot(1:max_tau,idx_d_s_sig,'r--');
    title(strcat(dataLabel,':Decision in the Past'));
    legend({'Pearson' 'Significance'});
    text(0.1,0.1,strcat('Max Tau:',num2str(idx_max_r_d_s)));
    ylim([0 1.1]);
    
    figure
    plot(1:max_tau,r_e_s,'b');
    hold on
    plot(1:max_tau,idx_e_s_sig,'r--');
    title(strcat(dataLabel,':Evidence in the Past'));
    legend({'Pearson' 'Significance'});
    text(0.1,0.1,strcat('Max Tau:',num2str(idx_max_r_e_s)));
    ylim([0 1.1]);
    
    tau_d = idx_max_r_d_s;
    
    figure
    plot(hEw((max_tau + 1):end),hDw((max_tau - tau_d + 1):(end - tau_d)),'b.');
    title(strcat(dataLabel,':Correlation'));
    xlabel('Evidence');
    ylabel('Decision - shifted');
 
    tau_e = idx_max_r_e_s;
    
    figure
    plot(hDw((max_tau + 1):end),hEw((max_tau - tau_e + 1):(end - tau_e)),'b.');
    title(strcat(dataLabel,':Correlation'));
    ylabel('Evidence');
    xlabel('Decision - shifted');
    
    
    
    function vector = setToColumnVector(vector)
        
        if size(vector,1) == 1
            
            vector = vector';
            
        end
        
    end

    
end

function [hEw, hDw] = setLocalMeanToZero(hEw,hDw)

    sizeWindow = 100;
    
    for i=1:sizeWindow:(length(hEw) - 2*sizeWindow)
            
        hEw(i:(sizeWindow + i - 1)) = zscore(hEw(i:(sizeWindow + i - 1)));
       
        hDw(i:(sizeWindow + i - 1)) = zscore(hDw(i:(sizeWindow + i - 1)));
  
    end
    
end

function plotPhaseDifferenceDistributionAtTransition(hEw,hDw,Trans,dataLabel)
        
        
    phases_e = hEw(Trans);
    phases_d = hDw(Trans);
    
    phases_diff = phases_e - phases_d;
    
    [unique_diff, nUnique_diff] = count_unique(phases_diff);
    
    figure
    bar(unique_diff,nUnique_diff);
    title(strcat(dataLabel,':Phases Differences at Transitions'));
        
end

function [hEw, hDw] = unwrapAllAlgorithm(hEw,hDw)
        
       
%         hEw = medfilt1(hEw);
%         hDw = medfilt1(hDw);
    
%         fhEw = fft(hEw);
%         fhDw = fft(hDw);
%         
%         L = length(hEw);
%         
%         H = [ones(1,(L/2)-1)*2 1 zeros(1,((L-1)-(L/2)+1))];
%         
%         ZfhEw = H.*fhEw;
%         ZfhDw = H.*fhDw;
%         
%         iZfhEw = ifft(ZfhEw);
%         iZfhDw = ifft(ZfhDw);
%         
%         hEw_chapeu = imag(iZfhEw);
%         hDw_chapeu = imag(iZfhDw);
%         
%         hE = real(iZfhEw);
%         hD = real(iZfhDw);

        ahE = hilbert(hEw);
        ahD = hilbert(hDw);
        
        rhE = real(ahE);
        rhD = real(ahD);
        
        ihE = imag(ahE);
        ihD = imag(ahD);
        
        hEw = applyKn(ihE,rhE);
        hDw = applyKn(ihD,rhD);
        
        function phiT = applyKn(phi_chapeu,phi)
 
            k(1) = 0;
            
            [pks,local_maxima] = findpeaks(phi_chapeu);
            
            [pks,local_minima] = findpeaks(-phi_chapeu);

%             dphi = diff(phi_chapeu);
%             idx_max = find(dphi>criterion);
%             idx_min = find(dphi<-criterion);
            
            for i=2:length(phi)
                
               if ~isempty(find(local_maxima == i))
                   
                   k(i) = k(i-1) + 1;
                   
               elseif ~isempty(find(local_minima == i))
                   
                   k(i) = k(i-1) - 1;
                   
               else
                   
                   k(i) = k(i-1);
                   
               end 
                
            end
           
            k = k';
            
            phiT = phi + 2*pi*k;
            
        end
        
end

function [hEw, hDw] = unwrapAll(hEw,hDw)

    hEw = unwrap(hEw);
    hDw = unwrap(hDw);

end

function [hEw, hDw] = smoothAll(hEw,hDw,span)

    hEw = smooth(hEw,span);

    hDw = smooth(hDw,span);

end

function getAndPlotConditionalPhaseDistribution(hEw,hDw,Trans,dataLabel)

    minTransInterval = min(diff(Trans));
    
    tau = ceil(linspace(1,minTransInterval,10));
    
    %for t=1:length(tau)
     for t=1:1   
       phase_e = hEw(Trans - tau(t)); 
       
       [unique_e, nUnique_e] = count_unique(phase_e);
       
       for j=1:length(unique_e)
           
           idx_e = find(hEw == unique_e(j));
           
           new_trans_idx = Trans - tau(t);
           
           idx_e = new_trans_idx(ismember(new_trans_idx,idx_e));
           
           phase_d = hDw(idx_e);
           
           [unique_d, nUnique_d] = count_unique(phase_d);
          
           figure
           bar(unique_d,nUnique_d);
           title(strcat(dataLabel));
           legend({strcat('Phase:',num2str(unique_e(j)))});
           
       end
        
    end



end

function plotFourierFrequenciesAmplitudes(hEw,hDw,dataLabel)

    fhEw = fft(hEw);
    fhDw = fft(hDw);
    
    
    afhEw = abs(fhEw);
    afhDw = abs(fhDw);
    
    afhEw(round(length(afhEw)/2):end) = [];
    afhDw(round(length(afhDw)/2:end)) = [];
    
    fq = (1:length(afhEw))./length(hEw);
    
    figure
    plot(fq,afhEw,'b');
    hold on
    plot(fq,afhDw,'r');
    legend({'Evidence' 'Decision'});
    title(strcat(dataLabel,':Fourier Spectrum'));
    
end

function plotOrbitsPhasesRelativeToTransitions(hEw,hDw,Trans,dataLabel)
        
    patch = -10:10;
    
%     ttE = hEw(patch + Trans(3)) - hDw(Trans(3));
%     ttD = hDw(patch + Trans(3)) - hDw(Trans(3));
%     
%     attE = round(radtodeg(ttE));
%     attD = round(radtodeg(ttD));
%     
%     dattE = diff(attE);
%     dattD = diff(attD);
%     
%     mdattE = max(dattE);
%     mdattD = max(dattD);
    
    figure
    for tt=1:length(Trans)

       plot(hEw(patch + Trans(tt)) - hDw(Trans(tt)),hDw(patch + Trans(tt)) - hDw(Trans(tt)),'b');
       hold on
       plot(hEw(Trans(tt)) - hDw(Trans(tt)),0,'ro');
       hold on

    end

    title(strcat(dataLabel,':Orbits - Phases Relative to Decision'));
    legend({'Phases' 'Transitions'});
    xlabel('Evidence');
    ylabel('Decision');
    
    figure
    for tt=1:length(Trans)

       plot(hDw(patch + Trans(tt)) - hEw(Trans(tt)),hEw(patch + Trans(tt)) - hEw(Trans(tt)),'b');
       hold on
       plot(hDw(Trans(tt)) - hEw(Trans(tt)),0,'ro');
       hold on

    end

    title(strcat(dataLabel,':Orbits - Phases Relative to Evidence'));
    legend({'Phases' 'Transitions'});
    ylabel('Evidence');
    xlabel('Decision');
        
end

function plotPhaseDistributionAtTau(hEw,hDw,Trans,dataLabel)
        
    Trans = Trans(10:end-10);
    
    mTransInterval = mean(diff(Trans));
    
    max_value = min((Trans(1)-1),(length(hEw)-Trans(end)-1));
    
    t = ceil(linspace(mTransInterval,max_value,10));
    
    i = 0;
    
    for tau=t
        
        i = i + 1;
    
        for tt=1:length(Trans)
            
            relPhaseE(tau - t(i) + i).phases(tt) = hEw(Trans(tt) + tau);
            relPhaseD(tau - t(i) + i).phases(tt) = hDw(Trans(tt) + tau);
        
        end
        
    end
    
    i = 0;
    
    for tau=t
       
        i = i + 1;
        
        for tt=1:length(Trans)
            
            phases_e(tt) = relPhaseE(tau - t(i) + i).phases(tt);
            phases_d(tt) = relPhaseD(tau - t(i) + i).phases(tt);

        end
        
        [uniques_e, numUnique_e] = count_unique(phases_e);
        [uniques_d, numUnique_d] = count_unique(phases_d);

        limit = max(uniques_e);
        
        figure
        bar(uniques_e,numUnique_e);
        xlim([-limit limit]);
        legend({'Evidence'});
        title(strcat('Phase Distribution at Tau: ',dataLabel));
        text(-limit,0,strcat('Tau:',int2str(t(i))));

        limit = max(uniques_d);
        
        figure
        bar(uniques_d,numUnique_d);
        xlim([-limit limit]);
        legend({'Decision'});
        title(strcat('Phase Distribution at Tau: ',dataLabel));
        text(-limit,0,strcat('Tau:',int2str(t(i))));
 
        
     end
    
end

function plotPhasesInTime(hEw,hDw,Trans,dataLabel)
        
        patch = -10:10;
        
        for tt=1:length(Trans)

            plot(patch,hEw(patch + Trans(tt)),'b');
            hold on
            plot(patch,hDw(patch + Trans(tt)),'r');
            hold on

        end
        
        title(strcat(dataLabel,':Phases'));
        legend({'Evidence' 'Decision'});
        xlabel('Time');
        ylabel('Phases');

end

function [rdtPCi, rdtPDi] = getDetrend(hEw,hDw,Trans)
        
    PCi = hEw;
    PDi = hDw;
    
    Tevent = Trans;
    
    fs = 14;
    
    % time axis

    Tstart = 1;
    Tend   = length(PCi);
    Ti     = 1:Tend;
    
    % number of transition events

    Nevent = length(Tevent);
    
    %% remove steps

    dPCi = diff( PCi );

    dPDi = diff( PDi );

%     % histogram of step sizes BEFORE removal
% 
%     ShowStepSizes( dPCi, dPDi );
%     title( 'before', 'FontSize', fs );

    % remove steps larger pi and smaller -pi

    kk = find( dPCi < -pi );
    dPCi(kk) = dPCi(kk) + 2*pi;

    jj = find( dPCi >  pi );
    dPCi(jj) = dPCi(jj) - 2*pi;

    NCstep = length(kk) + length(jj);

    kk = find( dPDi < -pi );
    dPDi(kk) = dPDi(kk) + 2*pi;

    jj = find( dPDi >  pi );
    dPDi(jj) = dPDi(jj) - 2*pi;

    NDstep = length(kk) + length(jj);

%     % histogram of step sizes AFTER removal
% 
%     ShowStepSizes( dPCi, dPDi );
%     title( 'after', 'FontSize', fs );

%     % compare number of removed jumps with number of events
% 
%     [Nevent NCstep NDstep]


    %% subtract average rate of phase increase (detrending)

    mdPCi = mean( dPCi );
    mdPDi = mean( dPDi );


    % reconstitute detrended signals

    rdtPCi = cumsum( [0 dPCi - mdPCi] );

    rdtPDi = cumsum( [0 dPDi - mdPDi] );
        
        
end

function getAndPlotDetrend(hEw,hDw,Trans,dataLabel)

    PCi = hEw;
    PDi = hDw;
    
    Tevent = Trans;
    
    fs = 14;
    
    % time axis

    Tstart = 1;
    Tend   = length(PCi);
    Ti     = 1:Tend;
    
    % number of transition events

    Nevent = length(Tevent);
    
    %% remove steps

    dPCi = diff( PCi );

    dPDi = diff( PDi );

%     % histogram of step sizes BEFORE removal
% 
%     ShowStepSizes( dPCi, dPDi );
%     title( 'before', 'FontSize', fs );

    % remove steps larger pi and smaller -pi

    kk = find( dPCi < -pi );
    dPCi(kk) = dPCi(kk) + 2*pi;

    jj = find( dPCi >  pi );
    dPCi(jj) = dPCi(jj) - 2*pi;

    NCstep = length(kk) + length(jj);

    kk = find( dPDi < -pi );
    dPDi(kk) = dPDi(kk) + 2*pi;

    jj = find( dPDi >  pi );
    dPDi(jj) = dPDi(jj) - 2*pi;

    NDstep = length(kk) + length(jj);

%     % histogram of step sizes AFTER removal
% 
%     ShowStepSizes( dPCi, dPDi );
%     title( 'after', 'FontSize', fs );

%     % compare number of removed jumps with number of events
% 
%     [Nevent NCstep NDstep]


    %% subtract average rate of phase increase (detrending)

    mdPCi = mean( dPCi );
    mdPDi = mean( dPDi );


    % reconstitute detrended signals

    rdtPCi = cumsum( [0 dPCi - mdPCi] );

    rdtPDi = cumsum( [0 dPDi - mdPDi] );


    % note that mean is not zero!  only final value is guaranteed to be zero
    % ...


    % reconstitute full signals

    rPCi = rdtPCi + Ti * mdPCi;

    rPDi = rdtPDi + Ti * mdPDi;


    figure;

    subplot(2,1,1);

    hold on;
    plot( Ti, rdtPCi, 'r', 'LineWidth', 1 );

    plot( Ti, rdtPDi, 'b', 'LineWidth', 1 );

    h = legend( 'C', 'D');
    set( h, 'FontSize', fs );

    plot( Ti(Tevent), rdtPCi(Tevent), 'r.' );
    plot( Ti(Tevent), rdtPDi(Tevent), 'b.' );

    hold off;

    xlabel('time steps', 'FontSize', fs );
    ylabel('detrended HT \Phi', 'FontSize', fs );

    subplot(2,1,2);

    hold on;
    plot( Ti, rPCi, 'r', 'LineWidth', 1 );

    plot( Ti, rPDi, 'b', 'LineWidth', 1 );

    h = legend( 'C', 'D');
    set( h, 'FontSize', fs );

    plot( Ti(Tevent), rPCi(Tevent), 'r.' );
    plot( Ti(Tevent), rPDi(Tevent), 'b.' );

    hold off;

    xlabel('time steps', 'FontSize', fs );
    ylabel('HT \Phi', 'FontSize', fs );
    
    title(dataLabel);


    %% plot orbits around selected transition events

    Nselect = 20;   % must be smaller than Nevent

    idx = randperm(Nevent);          

    Tselect = Tevent( idx(1:Nselect) );   % choose Nselect of Nevent events

    fullwindow = 51;
    halfwindow = 25;

    figure;
    
    subplot(1,2,1);

    hold on;
    for i=1:Nselect

        ievent = Tselect(i);

        iwindow = ievent+[-halfwindow:halfwindow];

        rdtPC = rdtPCi(iwindow) - rdtPCi(ievent);
        rdtPD = rdtPDi(iwindow) - rdtPDi(ievent);

        plot( rdtPC, rdtPD, 'r', 'LineWidth', 1 );
    end

    plot( 0, 0, 'k.' );

    hold off;

    axis 'square';

    xlabel('\Delta \Phi_C', 'FontSize', fs );
    ylabel('\Delta \Phi_D', 'FontSize', fs );

    title(strcat(dataLabel,':detrended'), 'FontSize', fs );

    subplot(1,2,2);

    hold on;
    for i=1:Nselect

        ievent = Tselect(i);

        iwindow = ievent+[-halfwindow:halfwindow];

        rPC = rPCi(iwindow) - rPCi(ievent);
        rPD = rPDi(iwindow) - rPDi(ievent);

        plot( rPC, rPD, 'r', 'LineWidth', 1 );
    end

    plot( 0, 0, 'k.' );

    hold off;

    axis 'equal';
    axis 'square';

    xlabel('\Delta \Phi_C', 'FontSize', fs );
    ylabel('\Delta \Phi_D', 'FontSize', fs );

    title(strcat(dataLabel,':reconstructed'), 'FontSize', fs );
    


end

function [hEw_fwht, hDw_fwht] = getAndPlotWalshHadamard(hEw,hDw,Trans,limit)
        
        hEw_fwht = fwht(hEw);
        hDw_fwht = fwht(hDw);

        hEw_fwht(limit+1:end) = [];
        hDw_fwht(limit+1:end) = [];

        figure
        plot(1:limit,hEw_fwht,'b');
        hold on
        plot(1:limit,hDw_fwht,'r');
        hold on
%         plot(Trans,ones(length(Trans),1),'ko');
        legend({'Evidence' 'Decision' 'Transitions'});
        title(strcat('Phases (',dataLabel,')'))       
        
end

function plotRelativePhaseDistributionAtTrantisionsWithTau(hEw,hDw,Trans,dataLabel)
        
       
    Trans = Trans(10:end-10);
    
    mTransInterval = mean(diff(Trans));
    
    max_value = min((Trans(1)-1),(length(hEw)-Trans(end)-1));
    
    t = ceil(linspace(mTransInterval,max_value,10));
    
    i = 0;
    
    for tau=t
        
        i = i + 1;
    
        for tt=1:length(Trans)
            
            relPhaseE(tau - t(i) + i,tt).phases = hEw - hEw(Trans(tt) + tau);
            relPhaseD(tau - t(i) + i,tt).phases = hDw - hDw(Trans(tt) + tau);
        
        end
        
    end
    
    i = 0;
    
    for tau=t
       
        i = i + 1;
        
        for tt=1:length(Trans)
            
            phases_e(tt) = relPhaseE(tau - t(i) + i,tt).phases(Trans(tt));
            phases_d(tt) = relPhaseD(tau - t(i) + i,tt).phases(Trans(tt));

        end
        
        [uniques_e, numUnique_e] = count_unique(phases_e);
        [uniques_d, numUnique_d] = count_unique(phases_d);

        limit = max(uniques_e);
        
        figure
        bar(uniques_e,numUnique_e);
        xlim([-limit limit]);
        legend({'Evidence'});
        title(strcat('Relative Phase Distribution at Transitions: ',dataLabel));
        text(-limit,0,strcat('Tau:',int2str(t(i))));

        limit = max(uniques_d);
        
        figure
        bar(uniques_d,numUnique_d);
        xlim([-limit limit]);
        legend({'Decision'});
        title(strcat('Relative Phase Distribution at Transitions: ',dataLabel));
        text(-limit,0,strcat('Tau:',int2str(t(i))));
 
        
    end
        
end

function getAndPlotOrbitsAndRelativePhases(hEw,hDw,Trans,dataLabel)

    Trans = Trans(10:end-10);
    
    mTransInterval = mean(diff(Trans));
    
    max_value = min((Trans(1)-1),(length(hEw)-Trans(end)-1));
    
    t = ceil(linspace(mTransInterval,max_value,10));
    
    i = 0;
    
    for tau=t
        
        i = i + 1;
    
        for tt=1:length(Trans)
            
            relPhaseE(tau - t(i) + i,tt).phases = hEw - hEw(Trans(tt) + tau);
            relPhaseD(tau - t(i) + i,tt).phases = hDw - hDw(Trans(tt) + tau);
        
        end
        
    end
    
    patch = -10:10;
    
    for i=1:length(t)
        
        figure
        for tt=1:length(Trans)

           plot(relPhaseE(i,tt).phases(patch + Trans(tt)),relPhaseD(i,tt).phases(patch + Trans(tt)),'b');
           hold on
           plot(relPhaseE(i,tt).phases(Trans(tt)),relPhaseD(i,tt).phases(Trans(tt)),'ro');
           hold on
 
           mRelPhaTransE(tt) = min(relPhaseE(i,tt).phases(patch + Trans(tt)));
           mRelPhaTransD(tt) = min(relPhaseD(i,tt).phases(patch + Trans(tt)));
           
        end
        
        mRelE = min(mRelPhaTransE);
        mRelD = min(mRelPhaTransD);
        
        title(strcat(dataLabel,':Orbits'));
        legend({'Relative Phases' 'Transitions'});
        xlabel('Evidence');
        ylabel('Decision');
        text(mRelE,mRelD,strcat('Tau:',int2str(t(i))));
 
    end
    
%     for i=1:length(t)
%         
%         figure
%         for tt=1:length(Trans)
% 
%             plot(patch,relPhaseE(i,tt).phases(patch + Trans(tt)),'b');
%             hold on
%             plot(patch,relPhaseD(i,tt).phases(patch + Trans(tt)),'r');
%             hold on
%             
%            mRelPhaTransE(tt) = min(relPhaseE(i,tt).phases(patch + Trans(tt)));
%            mRelPhaTransD(tt) = min(relPhaseD(i,tt).phases(patch + Trans(tt)));
% 
%         end
%         
%         mRelE = min(mRelPhaTransE);
%         mRelD = min(mRelPhaTransD);
%         
%         title(strcat(dataLabel,':relative Phases'));
%         legend({'Evidence' 'Decision'});
%         xlabel('Time');
%         ylabel('Phases');
%         text(mRelE,mRelD,strcat('Tau:',int2str(t(i))));
%  
%     end

end

function getAndPlotRelativePhases(hEw,hDw,Trans,dataLabel)

    Trans = Trans(10:end-10);
    
    mTransInterval = mean(diff(Trans));
    
    max_value = min((Trans(1)-1),(length(hEw)-Trans(end)-1));
    
    t = ceil(linspace(mTransInterval,max_value,10));
    
    i = 0;
    
    for tau=t
        
        i = i + 1;
    
        for tt=1:length(Trans)
            
            relPhaseE(tau - t(i) + i,tt).phases = hEw - hEw(Trans(tt) + tau);
            relPhaseD(tau - t(i) + i,tt).phases = hDw - hDw(Trans(tt) + tau);
        
        end
        
    end
    
    patch = -10:10;
    
    for i=1:length(t)
        
        figure
        for tt=1:length(Trans)

           plot(patch,relPhaseE(i,tt).phases(patch + Trans(i)),'b');
           hold on
           plot(0,relPhaseE(i,tt).phases(Trans(i)),'ro');
           hold on
 
           mRelPhaTrans(tt) = min(relPhaseE(i,tt).phases(patch + Trans(i)));
           
        end
        
        title(strcat(dataLabel,':Evidence'));
        legend({'Phases' 'Transitions'});
        text(min(patch),min(mRelPhaTrans),strcat('Tau:',int2str(t(i))));
        
        figure
        for tt=1:length(Trans)
             
           plot(patch,relPhaseD(i,tt).phases(patch + Trans(i)),'b');
           hold on
           plot(0,relPhaseD(i,tt).phases(Trans(i)),'ro');
           hold on
           
           mRelPhaTrans(tt) = min(relPhaseD(i,tt).phases(patch + Trans(i)));
        
        end
        
        title(strcat(dataLabel,':Decision'));
        legend({'Phases' 'Transitions'});
        text(min(patch),min(mRelPhaTrans),strcat('Tau:',int2str(t(i))));
 
    end

end

function getGranger(y,x,dataLabel)

    alpha = 0.05;
    max_lag = 5;

    [xw_yw_F, xw_yw_CV] = granger_cause(y,x,alpha,max_lag);
    [yw_xw_F, yw_xw_CV] = granger_cause(x,y,alpha,max_lag);

    %% Plots

    % Granger Causality between Hilbert Transform Phases with Transitions
    disp(strcat(dataLabel,':Decision Granger-cause Evidence:',int2str(xw_yw_F>xw_yw_CV)));
    disp(strcat('F-value:',num2str(xw_yw_F)));
    disp(strcat('Critical Value:',num2str(xw_yw_CV)));

    disp(strcat(dataLabel,':Evidence Granger-cause Decision:',int2str(yw_xw_F>yw_xw_CV)));
    disp(strcat('F-value:',num2str(yw_xw_F)));
    disp(strcat('Critical Value:',num2str(yw_xw_CV)));

end

function getAndPlotConditionalExpectations(hEw,hDw,Trans,dataLabel)
        
    meanInterval = round(2*mean(diff(Trans)));
    
    Trans(Trans<meanInterval) = [];
    Trans(Trans>(length(hEw)-meanInterval)) = [];
    
    max_tau = meanInterval;
    max_tau = 400;
    
    TEe_d = zeros(max_tau);
    TEd_e = zeros(max_tau);
    
    for tau_e = 1:max_tau
        
        for tau_d = tau_e:max_tau
        
            [idx_e_s_e, idx_e_i_e, idx_d_s_d, idx_d_i_d] = getIdxRelativePhases(hEw,hDw,Trans,tau_e,tau_d);
  
            TEe_d(tau_e,tau_d) = mean(hDw(Trans(idx_e_s_e) - tau_d)) - mean(hDw(Trans(idx_e_i_e) - tau_d));

            TEd_e(tau_e,tau_d) = mean(hEw(Trans(idx_d_s_d) - tau_e)) - mean(hEw(Trans(idx_d_i_d) - tau_e));
            
        end
    
    end
    
    TEe_e = zeros(1,max_tau);
    TEd_d = zeros(1,max_tau);
    
    for tau = 1:max_tau
        
            [idx_e_s_e, idx_e_i_e, idx_d_s_d, idx_d_i_d] = getIdxRelativePhases(hEw,hDw,Trans,tau,tau);
        
            TEe_e(tau) = mean(hEw(Trans(idx_e_s_e) - tau)) - mean(hEw(Trans(idx_e_i_e) - tau));

            TEd_d(tau) = mean(hDw(Trans(idx_d_s_d) - tau)) - mean(hDw(Trans(idx_d_i_d) - tau)); 
            
    end

    %% TE: Evidence (Phi) --> Evidence (Phi)
    figure
    plot(1:max_tau,TEe_e,'b.');
%     colorbar;
    title(strcat(dataLabel,': TE - \Phi\Phi'));
    xlabel('Tau');
    ylabel('Tau_\Phi in Time - Distance from Transition');
%     colormap(jet);
    ylim([min(TEe_e) max(TEe_e)]);
    
%     [x_zeros, y_zeros] = find(TEe_e == 0);
%     
%     x_zeros = x_zeros - size(TEe_e,1) + 1;
%     y_zeros = y_zeros - size(TEe_e,2) + 1;
%     
%     hold on     
%     for k = 1:numel(x_zeros)
%         
%         if y_zeros(k) > x_zeros(k)
%             
%             plot([x_zeros(k) x_zeros(k)],[y_zeros(k) y_zeros(k)], 'r.');
%             
%         end
%         
%     end
    
    %% TE: Decision (Psi) --> Decision (Psi)
    figure
    plot(1:max_tau,TEd_d,'b.');
%     colorbar;
    title(strcat(dataLabel,': TE - \Psi\Psi'));
    xlabel('Tau');
    ylabel('Tau_\Psi in Time - Distance from Transition');
%     colormap(jet);
    ylim([min(TEd_d) max(TEd_d)]);
    
%     [x_zeros, y_zeros] = find(TEd_d == 0);
%     
%     x_zeros = x_zeros - size(TEd_d,1) + 1;
%     y_zeros = y_zeros - size(TEd_d,2) + 1;
%     
%     hold on     
%     for k = 1:numel(x_zeros)
%         
%         if y_zeros(k) > x_zeros(k)
%             
%             plot([x_zeros(k) x_zeros(k)],[y_zeros(k) y_zeros(k)], 'r.');
%             
%         end
%         
%     end
    
    %% TE: Evidence (Phi) --> Decision (Psi)
    figure
    contour(TEe_d,min(min(TEe_d)):0.5:max(max(TEe_d)));
    colorbar;
    title(strcat(dataLabel,': TE - \Phi\Psi'));
    xlabel('Tau_\Phi in Time - Distance from Transition');
    ylabel('Tau_\Psi in Time - Distance from Transition');
    colormap(jet);
    
%     [x_zeros, y_zeros] = find(TEe_d == 0);
%     
%     x_zeros = x_zeros - size(TEe_d,1) + 1;
%     y_zeros = y_zeros - size(TEe_d,2) + 1;
%     
%     hold on     
%     for k = 1:numel(x_zeros)
%         
%         if y_zeros(k) > x_zeros(k)
%             
%             plot([x_zeros(k) x_zeros(k)],[y_zeros(k) y_zeros(k)], 'r.');
%             
%         end
%         
%     end
    
    %% TE: Decision (Psi) --> Evidence (Phi)
    figure
    contour(TEd_e,min(min(TEd_e)):0.5:max(max(TEd_e)));
    colorbar;
    title(strcat(dataLabel,': TE - \Psi\Phi'));
    xlabel('Tau_\Psi in Time - Distance from Transition');
    ylabel('Tau_\Phi in Time - Distance from Transition');
    colormap(jet);
    
%     [x_zeros, y_zeros] = find(TEd_e == 0);
%     
%     x_zeros = x_zeros - size(TEd_e,1) + 1;
%     y_zeros = y_zeros - size(TEd_e,2) + 1;
%     
%     hold on     
%     for k = 1:numel(x_zeros)
%         
%         if y_zeros(k) > x_zeros(k)
%             
%             plot([x_zeros(k) x_zeros(k)],[y_zeros(k) y_zeros(k)], 'r.');
%             
%         end
%         
%     end
    
%     tau = -max_tau:max_tau;
%     
%     figure
%     plot(tau,TEe_e,'b');
%     title(strcat(dataLabel,':Conditional Expectations'));
%     legend({'Evidence - Evidence'});
%     
%     figure
%     plot(tau,TEd_d,'b');
%     title(strcat(dataLabel,':Conditional Expectations'));
%     legend({'Decision - Decision'});
%         
%     figure
%     plot(tau,TEe_d,'b');
%     title(strcat(dataLabel,':Conditional Expectations'));
%     legend({'Evidence - Decision'});
%         
%     figure
%     plot(tau,TEd_e,'b');
%     title(strcat(dataLabel,':Conditional Expectations'));
%     legend({'Decision - Evidence'});
%     

    function [idx_e_s_e, idx_e_i_e, idx_d_s_d, idx_d_i_d] = getIdxRelativePhases(hEw,hDw,Trans,tau_e,tau_d)
        
            relPhaseE = hEw(Trans - tau_e) - hEw(Trans);
            relPhaseD = hDw(Trans - tau_d) - hDw(Trans);

            mrelPhaseE = mean(relPhaseE);

            mrelPhaseD = mean(relPhaseD);

            idx_e_s_e = find(relPhaseE >= mrelPhaseE);
            idx_e_i_e = find(relPhaseE < mrelPhaseE);

            idx_d_s_d = find(relPhaseD >= mrelPhaseD);
            idx_d_i_d = find(relPhaseD < mrelPhaseD);
            
    end

end

function [hEw, hDw] = absDegreePhases(hEw,hDw)
        
    % Transform Negative Phases in Positive Ones
    
    hEw = abs(hEw);
    hDw = abs(hDw);
    
    % Transform to Degrees (rounding)

    hEw = round(radtodeg(hEw));
    hDw = round(radtodeg(hDw));
    
end

function getAndPlotAutoAndCrossCorrelations(E,D,hEw,hDw,dataLabel)
    
    [e_auto, e_auto_lag] = xcorr(E);
    [ew_auto, ew_auto_lag] = xcorr(hEw);
    
    [d_auto, d_auto_lag] = xcorr(D);
    [dw_auto, dw_auto_lag] = xcorr(hDw);
    
    [ed_cross, ed_cross_lag] = xcorr(E,D);
    [edw_cross, edw_cross_lag] = xcorr(hEw,hDw);
    
    figure
    plot(e_auto_lag,e_auto,'b');
    title(strcat(dataLabel,':Auto-Correlation - Raw Data'));
    legend({'Evidence'});
    
    figure
    plot(d_auto_lag,d_auto,'r');
    title(strcat(dataLabel,':Auto-Correlation - Raw Data'));
    legend({'Decision'});
    
    figure
    plot(ew_auto_lag,ew_auto,'b');
    title(strcat(dataLabel,':Auto-Correlation - Phases'));
    legend({'Evidence'});
    
    figure
    plot(dw_auto_lag,dw_auto,'r');
    title(strcat(dataLabel,':Auto-Correlation - Phases'));
    legend({'Decision'});
    
    figure
    plot(ed_cross_lag,ed_cross,'b');
    title(strcat(dataLabel,':Cross-Correlation'));
    legend({'Evidence-Decision (raw data)'});
    
    figure
    plot(edw_cross_lag,edw_cross,'r');
    title(strcat(dataLabel,':Cross-Correlation'));
    legend({'Evidence - Decision (phases)'});
    
end

function [E_seg, D_seg] = getAndPlotSegmentsAroundTransitions(E,D,Trans,dataLabel)
        

    mean_half_int = mean(diff(Trans));
    
    Trans(Trans<mean_half_int) = [];
    Trans(Trans>(length(E)-mean_half_int)) = [];
    
    for i=1:length(Trans)
        
       E_seg(i).trans = Trans(i);
       D_seg(i).trans = Trans(i);
       
       E_seg(i).seg = E((Trans(i)-mean_half_int):(Trans(i)+mean_half_int));
       D_seg(i).seg = D((Trans(i)-mean_half_int):(Trans(i)+mean_half_int));
       
       [E_seg(i).hEa, E_seg(i).hEw, E_seg(i).hEwt, E_seg(i).hEy] = hilbertTransform(E_seg(i).seg);
       [D_seg(i).hDa, D_seg(i).hDw, D_seg(i).hDwt, D_seg(i).hDy] = hilbertTransform(D_seg(i).seg);
           
    end
    
    t = 3;
    
    figure
    plot(1:2*mean_half_int+1,E_seg(t).hEw,'b');
    hold on
    plot(1:2*mean_half_int+1,D_seg(t).hDw,'r');
    title(strcat(dataLabel,':Phase Segment'));
    legend({'Evidence' 'Decision'});
    
    figure
    plot(1:2*mean_half_int+1,E_seg(t).seg,'b');
    hold on
    plot(1:2*mean_half_int+1,D_seg(t).seg,'r');
    title(strcat(dataLabel,':Raw Data Segment'));
    legend({'Evidence' 'Decision'});
    
end

function plotPhaseDistributionAtTransitions(ew,dw,Trans,dataLabel)

    phases_e = ew(Trans);
    phases_d = dw(Trans);
    
    [uniques_e, numUnique_e] = count_unique(phases_e);
    [uniques_d, numUnique_d] = count_unique(phases_d);
    
    limit = max(uniques_e);
    
    figure
    bar(uniques_e,numUnique_e);
    xlim([-limit limit]);
    legend({'Evidence'});
    title(strcat('Phase Distribution at Transitions: ',dataLabel));
    
    limit = max(uniques_d);
    
    figure
    bar(uniques_d,numUnique_d);
    xlim([-limit limit]);
    legend({'Decision'});
    title(strcat('Phase Distribution at Transitions: ',dataLabel));
        
end

function plotHilbert(E,D,dr,er,dw,ew,dwt,ewt,Trans,dataLabel)

    figure
    plot(1:length(E),E,'b');
    hold on
    plot(1:length(D),D,'r');
    hold on
    plot(Trans,ones(length(Trans),1),'ko');
    legend({'Evidence' 'Decision'});
    title(strcat(dataLabel,':Orbits'));
    
%     figure
%     plot(1:vector_size,dr,'b');
%     hold on
%     plot(1:vector_size,dw,'r');
%     hold on
%     plot(1:vector_size,dwt,'g');
%     hold on
%     plot(Trans,ones(length(Trans),1),'ko');
%     legend({'Amplitude' 'Phase' 'Frequency' 'Transitions'});
%     title(strcat('Decision (',dataLabel,')'));
%     
%     figure
%     plot(1:vector_size,er,'b');
%     hold on
%     plot(1:vector_size,ew,'r');
%     hold on
%     plot(1:vector_size,ewt,'g');
%     hold on
%     plot(Trans,ones(length(Trans),1),'ko');
%     legend({'Amplitude' 'Phase' 'Frequency' 'Transitions'});
%     title(strcat('Evidence (',dataLabel,')'));
    
    figure
    plot(1:length(ew),ew,'b');
    hold on
    plot(1:length(dw),dw,'r');
    hold on
    plot(Trans,ones(length(Trans),1),'ko');
    legend({'Evidence' 'Decision' 'Transitions'});
    title(strcat('Phases (',dataLabel,')'))       
        
end

function [E, D, Trans, dataLabel] = getData(kind,dataLength)
        
    switch kind
    
        case 'Sanity' 
            
            %% Load DataSet Sanity
            data = load('sanity.mat');
            D = data.Di(1:dataLength);
            E = data.Ci(1:dataLength);
            Trans = data.Tevent(data.Tevent<dataLength);
            dataLabel = 'Sanity';
    
        case 'Robin'

            % Load DataSet Robin
            data = load('RobinTimeSeries.mat');
            trans = load('Time-Transitions.mat');
            D = data.X2dt(1:dataLength) - data.Y2dt(1:dataLength); %% Decision
            E = data.X1dt(1:dataLength) - data.Y1dt(1:dataLength); %% Evidence
            Trans = trans.trans.time(trans.trans.time<dataLength);
            dataLabel = 'Robins Model';
    
        case 'Random'

            D = rand(1,dataLength);
            E = rand(1,dataLength);
            transInterval = 2000;
            nTrans = ceil(dataLength/transInterval);
            Trans = randperm(dataLength,nTrans);
            Trans = Trans(Trans<(dataLength-transInterval));
            Trans = Trans(Trans>transInterval);
            dataLabel = 'Random';
            
         case 'Senoide'

            transInterval = 2000;
            nTrans = ceil(dataLength/transInterval);
            D = linspace(-nTrans*pi,nTrans*pi,transInterval*2*nTrans);
            E = sin(D);
            Trans = randperm(dataLength,nTrans);
            Trans = Trans(Trans<(dataLength-transInterval));
            Trans = Trans(Trans>transInterval);
            dataLabel = 'Senoide';
    
    end

end

function [amp, phase, freq, hW] = hilbertTransform(W)
        
       hW = hilbert(W);
     
       amp = abs(hW);
       phase = angle(hW);
       freq = [0 diff(phase)/2*pi]; 
        
end

function [uniques,numUnique] = count_unique(x,option)
%COUNT_UNIQUE  Determines unique values, and counts occurrences
%   [uniques,numUnique] = count_unique(x)
%
%   This function determines unique values of an array, and also counts the
%   number of instances of those values.
%
%   This uses the MATLAB builtin function accumarray, and is faster than
%   MATLAB's unique function for intermediate to large sizes of arrays for integer values.  
%   Unlike 'unique' it cannot be used to determine if rows are unique or 
%   operate on cell arrays.
%
%   If float values are passed, it uses MATLAB's logic builtin unique function to
%   determine unique values, and then to count instances.
%
%   Descriptions of Input Variables:
%   x:  Input vector or matrix, N-D.  Must be a type acceptable to
%       accumarray, numeric, logical, char, scalar, or cell array of
%       strings.
%   option: Acceptable values currently only 'float'.  If 'float' is
%           specified, the input x vector will be treated as containing
%           decimal values, regardless of whether it is a float array type.
%
%   Descriptions of Output Variables:
%   uniques:    sorted unique values
%   numUnique:  number of instances of each unique value
%
%   Example(s):
%   >> [uniques] = count_unique(largeArray);
%   >> [uniques,numUnique] = count_unique(largeArray);
%
%   See also: unique, accumarray

% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2009-03-17

testFloat = false;
if nargin == 2 && strcmpi(option,'float')
    testFloat = true;
end

nOut = nargout;
if testFloat
    if nOut < 2
        [uniques] = float_cell_unique(x,nOut);
    else
        [uniques,numUnique] = float_cell_unique(x,nOut);
    end
else
    try %this will fail if the array is float or cell
        if nOut < 2
            [uniques] = int_log_unique(x,nOut);
        else
            [uniques,numUnique] = int_log_unique(x,nOut);
        end
    catch %default to standard approach
        if nOut < 2
            [uniques] = float_cell_unique(x,nOut);
        else
            [uniques,numUnique] = float_cell_unique(x,nOut);
        end
    end
end

end

function [uniques,numUnique] = int_log_unique(x,nOut)
%First, determine the offset for negative values
minVal = min(x(:));

%Check to see if accumarray is appropriate for this function
maxIndex = max(x(:)) - minVal + 1;
if maxIndex / numel(x) > 1000
    error('Accumarray is inefficient for arrays when ind values are >> than the number of elements')
end

%Now, offset to get the index
index = x(:) - minVal + 1;

%Count the occurrences of each index value
numUnique = accumarray(index,1);

%Get the values which occur more than once
uniqueInd = (1:length(numUnique))';
uniques = uniqueInd(numUnique>0) + minVal - 1;

if nOut == 2
    %Trim the numUnique array
    numUnique = numUnique(numUnique>0);
end
end 

function [uniques,numUnique] = float_cell_unique(x,nOut)

if ~iscell(x)
    %First, sort the input vector
    x = sort(x(:));
    numelX = numel(x);
    
    %Check to see if the array type needs to be converted to double
    currClass = class(x);
    isdouble = strcmp(currClass,'double');
    
    if ~isdouble
        x = double(x);
    end
    
    %Check to see if there are any NaNs or Infs, sort returns these either at
    %the beginning or end of an array
    if isnan(x(1)) || isinf(x(1)) || isnan(x(numelX)) || isinf(x(numelX))
        %Check to see if the array contains nans or infs
        xnan = isnan(x);
        xinf = isinf(x);
        testRep = xnan | xinf;
        
        %Remove all of these from the array
        x = x(~testRep);
    end
    
    %Determine break locations of unique values
    uniqueLocs = [true;diff(x) ~= 0];
else
    isdouble = true; %just to avoid conversion on finish
    
    %Sort the rows of the cell array
    x = sort(x(:));
    
    %Determine unique location values
    uniqueLocs = [true;~strcmp(x(1:end-1),x(2:end)) ~= 0] ;
end

%Determine the unique values
uniques = x(uniqueLocs);

if ~isdouble
    x = feval(currClass,x);
end

%Count the number of duplicate values
if nOut == 2
    numUnique = diff([find(uniqueLocs);length(x)+1]);
end


end


end

