
function h = binOptimalPSTHkernelGratings(condition,stimIds,spike_times,start_time,end_time)
  
   trials_label = find(stimIds == condition); 
   
   spike_times = spike_times./32000;
    
   trials_spikes = spike_times(trials_label,:);

   nTrials = size(trials_spikes,1);

   spikes_vector = reshape(trials_spikes.',[],1);

   spikes_vector = spikes_vector.';

   spikes_vector = sort(spikes_vector);

   spikes_vector = spikes_vector(spikes_vector>0);

   spikes_vector = spikes_vector(spikes_vector>=(start_time/1000) & spikes_vector<=(end_time/1000));

   NbinMax = round((end_time/1000 - start_time/1000) * 1000/2);

   bandwidths = (2/1000):(2/1000):(end_time/1000/2);

   [optN, psth.C, psth.N] = sshist(spikes_vector,2:NbinMax);

   [optW, kernel.C, kernel.W] = sskernel(spikes_vector,bandwidths);

   histData.spikes_vector = spikes_vector;
   histData.psth.optimalNBins = optN;
   histData.psth.Cost = psth.C;
   histData.psth.NBinsTested = psth.N;

   histData.kernel.optimalBinWidth = optW;
   histData.kernel.Cost = kernel.C;
   histData.kernel.WBinsTested = kernel.W;

   histData.nTrials = nTrials;

   binSize = (end_time/1000 - start_time/1000)/optN;

   for k=1:optN

       spikes = length(spikes_vector(spikes_vector>=((k-1)*binSize + start_time/1000) & spikes_vector<(k*binSize + start_time/1000)));

       binHeight(k) = spikes/(nTrials*binSize);

   end

   histData.psth.binHeight = binHeight;
   histData.psth.binSize = binSize;

   histData.psth.binHeightNormal = histData.psth.binHeight./max(histData.psth.binHeight);
    
   [histData.kernel.KY,histData.kernel.KX] = ksdensity(histData.spikes_vector,'width',histData.kernel.optimalBinWidth);
   
   more = length(histData.kernel.KX(histData.kernel.KX>end_time/1000));
   less = length(histData.kernel.KX(histData.kernel.KX<start_time/1000));
   histData.kernel.KX = histData.kernel.KX(less+1:end-more);
   histData.kernel.KY = histData.kernel.KY(less+1:end-more);
   
   f = figure;
   subplot(2,1,1);
   
   for k=1:size(trials_spikes,1)
       
       spikes = trials_spikes(k,:);
       spikes = spikes(spikes>0);
       spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
       spikes = sort(spikes);
       
       for i=1:size(spikes,2)

           plot([spikes(i) spikes(i)],[(k-0.4) (k)],'b');
           hold on;

       end
   
   end
   
   subplot(2,1,2);

   plot(histData.kernel.KX,histData.kernel.KY);
   
   g = figure;
   
   subplot(2,1,1);
   
   for i=1:size(spikes_vector,2)

        plot([spikes_vector(i) spikes_vector(i)],[(1-0.4) (1)],'b');
        hold on;

   end
       
   subplot(2,1,2);

   plot(histData.kernel.KX,histData.kernel.KY);
   
%    NBins = histData.psth.optimalNBins;
% 
%    KX = KX.*(NBins/5);
%    KY = KY./max(KY);
% 
%    f = figure;
%    [AX, H1, H2] = plotyy(1:NBins,histData.psth.binHeightNormal,KX,KY,@bar,@plot);
% 
%    set(AX,'XLim',[1 length(histData.psth.binHeightNormal)]);
%    set(H2,'Color','r');

   print(f,'-dbmp',strcat('grating-raster-all-condition-',int2str(condition),'.bmp'));
   print(g,'-dbmp',strcat('grating-raster-condensed-condition-',int2str(condition),'.bmp'));
   
   save(strcat('grating-condition-',int2str(condition)),'histData');
    
end