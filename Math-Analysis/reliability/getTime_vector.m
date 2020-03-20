function Time_Vectors = getTime_vector(trials_spikes,nTrials,HSIZE,sigma,start_time,end_time,Filter_Dimension,latency)
    
%%%   BEGIN TRIALS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:nTrials
       
       spikes = trials_spikes(i,:);
       
       spikes = spikes - latency;

       spikes = spikes(spikes>0);

       spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

       trials(i).spikes = sort(spikes);
       
      if Filter_Dimension == 0
 
          trials(i).spikes_conv = trials(i).spikes;

      elseif Filter_Dimension == 1
           
%            x = linspace(-HSIZE / 2, HSIZE / 2, HSIZE);
%            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
           
           hgauss = fspecial('gaussian',[HSIZE 1],sigma);
           
           trials(i).spikes_conv = imfilter(trials(i).spikes,hgauss);
           
       elseif Filter_Dimension == 2
           
           hgauss = fspecial('gaussian',[HSIZE HSIZE],sigma);
           
           trials(i).spikes_conv = imfilter(trials(i).spikes,hgauss);
           
       end
       
       expoente = 4;
       
       trials(i).spikesmiliseconds = round(trials(i).spikes_conv*10^expoente);
       
       trials(i).spikesmiliseconds = sort(trials(i).spikesmiliseconds);
       
    end

    
    for i=1:nTrials

        if length(trials(i).spikesmiliseconds) == 0
            
            sizes(i) = 0;
            
        else
            
            sizes(i) = max(trials(i).spikesmiliseconds);
            
        end
        
    end
    
    maximum = max(sizes);
    
    for i=1:nTrials
        
       trials(i).time_vector = zeros(1,maximum);

       nSpikes = length(trials(i).spikesmiliseconds);
       
       if nSpikes > 0
           
           for k=1:nSpikes

               if int64(trials(i).spikesmiliseconds(k)) > 0
                   
                    trials(i).time_vector(int64(trials(i).spikesmiliseconds(k))) = int64(trials(i).spikesmiliseconds(k));
                    
               end
               
           end
       
       end
      
       
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Vectors = trials;

end

