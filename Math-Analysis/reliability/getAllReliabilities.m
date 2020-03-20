function [TV, real] = getAllReliabilities(HSIZE,bandwidths,trials_spikes,nTrials,start_time,end_time,Filter_Dimension,how_to_do)

Time_Vectors = zeros(1,nTrials);
TV(1:length(bandwidths)) = struct('Time_Vectors',Time_Vectors);

parfor w=1:length(bandwidths)
        
    sigma = bandwidths(w);
    
    if strcmp(how_to_do,'one_by_one')

        HSIZE_ = HSIZE;

    elseif strcmp(how_to_do,'all_together')
        
        HSIZE_ = HSIZE(w);

    end

%   [a, b] = doTheThing(trials_spikes,nTrials,HSIZE_,sigma,start_time,end_time,Filter_Dimension);
%   Time_Vectors(w).trials = a;
%   real(w) = b;

    Time_Vector = getTime_vector(trials_spikes,nTrials,HSIZE_,sigma,start_time,end_time,Filter_Dimension);
    
    real(w) = getReliability(Time_Vector,nTrials);
   
    TV(w).Time_Vectors = Time_Vector;
        
end
   


% function [Time_Vector, reliability] = doTheThing(trials_spikes,nTrials,HSIZE_,sigma,start_time,end_time,Filter_Dimension)
%         
%     Time_Vector = getTime_vector(trials_spikes,nTrials,HSIZE_,sigma,start_time,end_time,Filter_Dimension);
% 
%     reliability = getReliability(Time_Vector,nTrials);
%     
% end
    
    
   
end

