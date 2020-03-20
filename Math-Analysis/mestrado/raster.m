function h = raster(spike_times)

    spike_times = spike_times./ 32000;
    l = size(spike_times,1);
    
    for i=1:l
        
        spike_times(i,:) = spike_times(i,spike_times(i,:)>0);
        
    end
    
    for i=1:l
        
        g = size(spike_times,2);
        
        for i=1:g

            plot([spike_times(l,g) list(l,g)],[(idx-0.8) (idx)],'r');
            hold on;

        end

        idx = idx + 1;
  
    end
  

end