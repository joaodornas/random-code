function h_tf_win_con_sanity

san = load('sanity.mat');

dataLength = 100000;

Ci = san.Ci(1:dataLength);
Di = san.Di(1:dataLength);
Trans = san.Tevent;

gc_windowing(Ci,Di,{'Ci','Di'},Trans);
%gc_windowing(Di,Ci,{'Di','Ci'},Trans);

function gc_windowing(T,S,label,trans)

    G = T./S;
    
    G(isinf(G)) = 0;
    G(isnan(G)) = 0;
    
    TS = fft(G);
    
    realTS = real(TS);
    analyticTS = hilbert(imag(TS));
    nhiTS = 1.*imag(analyticTS);
    
%     figure
%     plot(1:length(T),imTS,'b');
%     legend({'Imaginary'});
%     title(strcat('Target:',label{1}));
%     %hold on
%     figure
%     plot(1:length(T),nhrTS,'r');
%     legend({'Hilbert-Real'});
%     title(strcat('Target:',label{1}));
    
%    residuals = imTS - nhrTS;

    residuals = realTS - nhiTS;
    
    h = lillietest(residuals);
    
    disp(strcat('Is it Normally Distributed?:',int2str(1-h)));
    
    if h == 0
        
        [H, p] = ttest(residuals);
        
        disp('parametric');
        
    elseif h == 1
        
        [p, H] = signrank(residuals);
        
        disp('non-parametric');
        
    end

    disp(strcat('Is there causallity?:',int2str(1-H)));
    disp(strcat('p-value:',num2str(p)));
    
    disp(strcat('Mean:',num2str(mean(residuals))));
    disp(strcat('Zeros:',int2str(length(residuals(residuals==0)))));
    
    
%     % number of repetitions of analysis
%     n_repeats = 100;
% 
%     % number of selected transitions for each repetition, equals number of concatenated samples 
%     n_selected_transitions = 20;
% 
%     % number and index of samples in a window
%     s_range = 60;
%     s_spacing = 2;
%     ix_samples = 0:s_spacing:s_range;
%     n_samples = length(ix_samples);
% 
%     % number and index of offsets (at which samples are taken and evaluated)
%     o_range = 200;
%     o_spacing = 5;
%     ix_offset = -o_range:o_spacing:o_range;
%     n_offsets = length(ix_offset);
% 
%     % granger params
%     alpha = 0.005;
%     max_lag = 5;
%     
%     % center of transitions
%     transitions = trans;
%     
%     % loop over repetitions
%     for ir=1:n_repeats
% 
%         if mod(ir,20)==0
%             ir
%         end
% 
%         % number of transitions
%         nTrans = length(transitions);
%         
%         rand_trans = transitions(randperm(nTrans));   % only one argument to randperm!
% 
%         few_trans = rand_trans(1:n_selected_transitions);
% 
%         % indices to samples from all selected transitions (for offset zero)
%         % length is n_samples * n_selected_transitions;
%         
%         ix_as = sort( repmat(few_trans,1,n_samples) );
% 
%         ix_as = ix_as + repmat( ix_samples, 1, n_selected_transitions );
% 
%         % loop over offsets
%         for io=1:n_offsets
% 
%             T_w = T(ix_as + ix_offset(io));
%             S_w = S(ix_as + ix_offset(io));
%           
%             T_w = T_w - mean(T_w); 
%             S_w = S_w - mean(S_w); 
% 
%             T_w = T_w/std(T_w);
%             S_w = S_w/std(S_w);
% 
%             [GC_S_T.F(ir,io), GC_S_T.c_v(ir,io)] = granger_cause(T_w,S_w,alpha,max_lag);
% 
%         end
% 
%     end
% 
%     GC_S_T.F(isnan(GC_S_T.F)) = 0;
%     
%     mean_GC_S_T_F = mean(GC_S_T.F,1);
%  
%     f = figure;
%     plot(ix_offset,mean_GC_S_T_F,'r');
%     legend(label{2});
%     title(strcat('Target:',label{1}));

end

end

