
Folder.Settings= 'Settings/'; % path to settings folder 
Folder.Results= 'Results/'; % path to results folder 
Folder.Events= 'Events/'; % path to events folder 

SettingsFile = input('Please give the Settings file name : ', 's');

Settings= CExperimentalSettings(SettingsFile, Folder.Settings);

Results = load(strcat(Folder.Results,SettingsFile));
Blocks = Results.Block;
nBlocks = length(Blocks);

ExperimentResponse = zeros(Settings.Trials,nBlocks,3);

FrameRate = Settings.Screen.ScreenMode.hz;

offset = 5;

ihs = 1;
for j=1:nBlocks
   
    Sequence_Of_Trials = Blocks{j}.Sequence_Of_Trials;
    
    for i=1:length(Sequence_Of_Trials)
        
        if ~isempty(Blocks{j}.Trial(i).LostTrackFrame); 
            if Blocks{j}.Trial(i).LostTrackFrame <= (FrameRate/2) 
                hs = 1; 
            else
                hs = floor(Blocks{j}.Trial(i).LostTrackFrame/(FrameRate/2)) + 1;
            end
        else
            hs = 0;
        end
            
        ExperimentResponse(Sequence_Of_Trials(i),j,1) = hs;
        if isempty(Blocks{j}.Trial(i).LostTrackWasCorrect)
            LostTrackWasCorrect = 0; 
        else
            LostTrackWasCorrect = Blocks{j}.Trial(i).LostTrackWasCorrect; 
        end
        ExperimentResponse(Sequence_Of_Trials(i),j,2) = LostTrackWasCorrect;
        ExperimentResponse(Sequence_Of_Trials(i),j,3) = Blocks{j}.ProbeIsTarget(i);
        
        allHS(ihs) = hs;
        
        ihs = ihs + 1;
        
    end

end

figure

for i=1:length(Sequence_Of_Trials)
   
    subplot(5,1,i)
    
    events = load(strcat(Folder.Events,'eventsBalls-',SettingsFile,'-T',int2str(i)));
    
    eventsArea = events.eventsArea;
    idx_events_Area = find(eventsArea);
    idx_events_Area(idx_events_Area > ( max(allHS) + offset)) = [];
    
    eventsBallColorMatch = sum(events.eventsBallColorMatch,1)./3;
    idx_events_BallColorMatch = find(eventsBallColorMatch);
    idx_events_BallColorMatch(idx_events_BallColorMatch > ( max(allHS) + offset)) = [];
    
    eventsBallDistance = sum(events.eventsBallDistance,1)./3;
    idx_events_BallDistance = find(eventsBallDistance);
    idx_events_BallDistance(idx_events_BallDistance > ( max(allHS) + offset)) = [];
    
    eventsBallMixing = squeeze(sum(squeeze(sum(events.eventsBallMixing(:,:,1:end),1)),1))./3;
    idx_events_BallMixing = find(eventsBallMixing);
    idx_events_BallMixing(idx_events_BallMixing > ( max(allHS) + offset)) = [];
    
    plot(idx_events_Area,eventsArea(idx_events_Area),'r*');
    hold on
    plot(idx_events_BallColorMatch,eventsBallColorMatch(idx_events_BallColorMatch) + 1,'b*');
    hold on
    plot(idx_events_BallDistance,eventsBallDistance(idx_events_BallDistance) + 2,'g*');
    hold on
    plot(idx_events_BallMixing,eventsBallMixing(idx_events_BallMixing) + 3,'m*');
    hold on
    for j=1:nBlocks
        if ~ExperimentResponse(i,j,2) && (ExperimentResponse(i,j,3) == 1) 
            hs = ExperimentResponse(i,j,1);
            plot([hs hs],[0 6],'k');
            hold on
            plot([hs hs],[6 6],'ko');
        end
    end
    
    ylim([0 7]);
    title({strcat('Trial:',int2str(i))});
    
end

legend({'Area' 'Color Match' 'Displacement' 'Mixing' 'Lose Track'});

