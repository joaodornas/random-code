function reliability = getReliability(all_trials)

if islogical(all_trials)
    
    all_trials = double(all_trials);
    
end

%%%   CALCULATE RELIABILITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % disp('Calculate Reliability...');
    
    reliability = 0;
    
    totalTrials = size(all_trials,1);
    
    for i=1:totalTrials
       
        for j=(i+1):totalTrials
        
            normai = norm(all_trials(i,:));
            normaj = norm(all_trials(j,:));
            
            reliability = reliability + dot(all_trials(i,:),all_trials(j,:))/(normai * normaj);
            
        end
        
    end

    reliability = reliability .* 2 ./ (totalTrials.*(totalTrials - 1));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end