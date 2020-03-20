
nROI = 90;
nRuns = 4*6;

ireal = 0;

for iROI=1:nROI
   
    %allruns = BothTrackPassive(iROI,:);
    %allruns = TrackPassiveThre20(iROI,:);
    %allruns = OnlyTrack15(iROI,:);
    allruns = BothTrackPassive(iROI,:);
    
%     if sum(allruns) == nRuns
%         
%         ireal = ireal + 1;
%         
%         reliability(ireal) = iROI;
%         
%     end

    reliability(iROI) = sum(allruns)/nRuns;
    
end

reliability = reliability';