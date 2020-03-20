

disp('Loading Data');

nRuns = 32;

all_Track_pos = [];
all_Track_neg = [];
for iRun=1:nRuns
     load(strcat('Global-Track-Run-',int2str(iRun),'.mat'));
     all_Track_pos = [all_Track_pos; Track.run_pos];
     all_Track_neg = [all_Track_neg; Track.run_neg];
     clear Track
end

all_Passive_pos = [];
all_Passive_neg = [];
for iRun=(nRuns+1):nRuns*2
     load(strcat('Global-Passive-Run-',int2str(iRun),'.mat'));
     all_Passive_pos = [all_Passive_pos; Passive.run_pos];
     all_Passive_neg = [all_Passive_neg; Passive.run_neg];
     clear Passive
end

all_Rest_pos = [];
all_Rest_neg = [];
for iRun=(nRuns*2+1):nRuns*3
     load(strcat('Global-Rest-Run-',int2str(iRun),'.mat'));
     all_Rest_pos = [all_Rest_pos; Rest.run_pos];
     all_Rest_neg = [all_Rest_neg; Rest.run_neg];
     clear Rest
end

nTotalVoxels = 160990;
nROI = 90;
