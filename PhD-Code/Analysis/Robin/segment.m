function W = segment(data)

Robin = load(data);

    X(1,:) = Robin.X1dt(1:end);
    X(2,:) = Robin.Y1dt(1:end);
    X(3,:) = Robin.X2dt(1:end);
    X(4,:) = Robin.Y2dt(1:end);
    
Decision = X(3,:,1);
 
time = transitions(Decision);
 
w = 1;

for i=1:(length(time) -1)
    
   present_transition = time(i);
   next_transition = time(i+1);
   
   interval_between_transitions(i) = next_transition - present_transition;
   
end

mean_interval_time = floor(mean(interval_between_transitions));

three_period = mean_interval_time*3;

for i=1:length(time)
   
   if ( (time(i) - three_period) >= 1 ) && ( (time(i) + three_period) <= time(end) )
    
        past_EX(w,:) = X(1,(time(i) - three_period):time(i));
        past_EY(w,:) = X(2,(time(i) - three_period):time(i));
        past_DX(w,:) = X(3,(time(i) - three_period):time(i));
        past_DY(w,:) = X(4,(time(i) - three_period):time(i));
        
        transition_index = intersect((time(i) - three_period):time(i),time);
        
        this_period_vector = (time(i) - three_period):time(i);
        
        for j=1:length(transition_index)
            
            W.transitions_past(w).new_scale(j) = find(this_period_vector==transition_index(j));
            
        end
  
        future_EX(w,:) = X(1,(time(i):(time(i) + three_period)));
        future_EY(w,:) = X(2,(time(i):(time(i) + three_period)));
        future_DX(w,:) = X(3,(time(i):(time(i) + three_period)));
        future_DY(w,:) = X(4,(time(i):(time(i) + three_period)));
        
        transition_index = intersect(time(i):(time(i) + three_period),time);
        
        this_period_vector = time(i):(time(i) + three_period);
        
        for j=1:length(transition_index)
            
            W.transitions_future(w).new_scale(j) = find(this_period_vector==transition_index(j));
            
        end
       
        w = w + 1;
   
   end
  
end

    W.past_EX = past_EX;
    W.past_EY = past_EY;
    W.past_DX = past_DX;
    W.past_DY = past_DY;
    
    W.mean_past_EX = mean(past_EX,1);
    W.mean_past_EY = mean(past_EY,1);
    W.mean_past_DX = mean(past_DX,1);
    W.mean_past_DY = mean(past_DY,1);
    
    W.future_EX = future_EX;
    W.future_EY = future_EY;
    W.future_DX = future_DX;
    W.future_DY = future_DY;
    
    W.mean_future_EX = mean(future_EX,1);
    W.mean_future_EY = mean(future_EY,1);
    W.mean_future_DX = mean(future_DX,1);
    W.mean_future_DY = mean(future_DY,1);

    W.n_transitions = length(time);

    W.mean_interval_time = mean_interval_time;
    
    save('segment.mat','W');
    
end

