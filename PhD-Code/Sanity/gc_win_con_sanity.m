function gc_win_con_sanity

san = load('sanity.mat');

Ci = san.Ci;
Di = san.Di;
Trans = san.Tevent;

% gc_windowing(Ci,Di,{'Ci','Di'},Trans);
% gc_windowing(Di,Ci,{'Di','Ci'},Trans);

hCi = hilbert(Ci);
hDi = hilbert(Di);

amp_Ci = abs(hCi);
amp_Di = abs(hDi);

phase_Ci = angle(hCi);
phase_Di = angle(hDi);

gc_windowing(amp_Ci,amp_Di,{'Amp_Ci','Amp_Di'},Trans);
gc_windowing(amp_Di,amp_Ci,{'Amp_Di','Amp_Ci'},Trans);

gc_windowing(phase_Ci,phase_Di,{'Phase_Ci','Phase_Di'},Trans);
gc_windowing(phase_Di,phase_Ci,{'Phase_Di','Phase_Ci'},Trans);


function gc_windowing(T,S,label,trans)

    % number of repetitions of analysis
    n_repeats = 100;

    % number of selected transitions for each repetition, equals number of concatenated samples 
    n_selected_transitions = 20;

    % number and index of samples in a window
    s_range = 60;
    s_spacing = 2;
    ix_samples = 0:s_spacing:s_range;
    n_samples = length(ix_samples);

    % number and index of offsets (at which samples are taken and evaluated)
    o_range = 200;
    o_spacing = 5;
    ix_offset = -o_range:o_spacing:o_range;
    n_offsets = length(ix_offset);

    % granger params
    alpha = 0.005;
    max_lag = 5;
    
    % center of transitions
    transitions = trans;
    
    % loop over repetitions
    for ir=1:n_repeats

        if mod(ir,20)==0
            ir
        end

        % number of transitions
        nTrans = length(transitions);
        
        rand_trans = transitions(randperm(nTrans));   % only one argument to randperm!

        few_trans = rand_trans(1:n_selected_transitions);

        % indices to samples from all selected transitions (for offset zero)
        % length is n_samples * n_selected_transitions;
        
        ix_as = sort( repmat(few_trans,1,n_samples) );

        ix_as = ix_as + repmat( ix_samples, 1, n_selected_transitions );

        % loop over offsets
        for io=1:n_offsets

            T_w = T(ix_as + ix_offset(io));
            S_w = S(ix_as + ix_offset(io));
          
            T_w = T_w - mean(T_w); 
            S_w = S_w - mean(S_w); 

            T_w = T_w/std(T_w);
            S_w = S_w/std(S_w);

            [GC_S_T.F(ir,io), GC_S_T.c_v(ir,io)] = granger_cause(T_w,S_w,alpha,max_lag);

        end

    end

    GC_S_T.F(isnan(GC_S_T.F)) = 0;
    
    mean_GC_S_T_F = mean(GC_S_T.F,1);
 
    f = figure;
    plot(ix_offset,mean_GC_S_T_F,'r');
    legend(label{2});
    title(strcat('Target:',label{1}));

end

end

