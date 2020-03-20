function [m, A, preferida, anti_preferida, rate, DI, response, sd] = getDirectionalSelectivity(DTC,start_time,end_time)


Spass = load(DTC);

nConditions = 16;

angles = 22.5;

bin_size = (end_time - start_time) / 1000;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

k = 1;

for i=1:nConditions
    
    trials_spikes = spike_times(trials_label(i).labels(:),:);
    
    nTrials = size(trials_spikes,1);
    
    for n=1:nTrials
        
       spike_train = trials_spikes(n,:);
       
       spike_train = spike_train(spike_train>0);
        
       A_train = spike_train(spike_train>=start_time/1000 & spike_train<=end_time/1000);
       
       all_A(k) = numel(A_train)/( bin_size );
       
       ae_train = spike_train(spike_train>=0/1000 & spike_train<=start_time/1000);
       
       a_e(n) = numel(ae_train)/( bin_size );
       
       trial_train = spike_train(spike_train>=start_time/1000 & spike_train<=end_time/1000);
       
       trial(n) = numel(trial_train)/( bin_size );
       
       k = k + 1;
        
    end
    
     sd(i) = std(trial);
    
     cond(i).ae = mean(a_e);
    
     cond(i).rate = mean(trial);
    
     spike_train = reshape(trials_spikes.',[],1);
        
     spike_train = spike_train(spike_train>0);
     
     ae = spike_train(spike_train>=0/1000 & spike_train<=start_time/1000);
        
     spike_train = spike_train(spike_train>=start_time/1000 & spike_train<=end_time/1000);
     
     all_m(i) = numel(ae)/( (nTrials) * bin_size );
     
     rate(i) = numel(spike_train)/( (nTrials) * bin_size );
     
     clear trial;
    
end

all_ae = 0;

for i=1:nConditions
    
    all_ae = all_ae + cond(i).ae;
    
end

all_ae = all_ae / nConditions;

idx = find(rate==max(rate));

preferida = (idx - 1) * angles;

anti_preferida = preferida + 180;

if anti_preferida > 360
    
    anti_preferida = anti_preferida - 360;
    
end

anti_preferida_idx = ( anti_preferida / angles ) + 2;

if anti_preferida_idx > nConditions, anti_preferida_idx = anti_preferida_idx- nConditions; end

A = max(all_A);

m = mean(all_m);


DI = getDI(rate(anti_preferida_idx),m,rate(idx));

DI_trials = getDI(cond(anti_preferida_idx).rate,all_ae,cond(idx).rate);

if DI >= 0.5
    
    response = 'direction selective';
    
elseif DI < 0.5
    
    response = 'non-direction selective';
    
end

function DI = getDI(R_anti, R_spont, R_pref)

    
    DI = 1 - ( ( R_anti - R_spont ) / ( R_pref - R_spont ) ) ;
    

end

end

