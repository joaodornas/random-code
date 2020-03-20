
nRuns = 4;
nSubjects = 8;
nTotalRuns = 32;

SUBJECT_FOLDER{1} = 'SUBJECT-1-22-10-2015';
SUBJECT_FOLDER{2} = 'SUBJECT-2-26-10-2015';
SUBJECT_FOLDER{3} = 'SUBJECT-3-3-11-2015';
SUBJECT_FOLDER{4} = 'SUBJECT-4-2-11-2015';
SUBJECT_FOLDER{5} = 'SUBJECT-5-2-11-2015';
SUBJECT_FOLDER{6} = 'SUBJECT-6-24-11-2015';
SUBJECT_FOLDER{7} = 'SUBJECT-7-14-01-2016';
SUBJECT_FOLDER{8} = 'SUBJECT-8-14-01-2016';

areaLabel{1} = 'Angular-L';
areaLabel{2} = 'Angular-R';
areaLabel{3} = 'Parietal-Sup-L';
areaLabel{4} = 'Parietal-Sup-R';
areaLabel{5} = 'Parietal-Inf-L';
areaLabel{6} = 'Parietal-Inf-R';

for iArea=1:length(areaLabel)
    
    disp(areaLabel{iArea});
    
    this_area = areaLabel{iArea};
    
    %% RESTING STATE
    
%     iiRun = 0;
%     for iSubject=1:nSubjects
% 
%         disp(int2str(iSubject));
%         
%         subject = SUBJECT_FOLDER{iSubject};
%         
%         load(strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\',subject,'\output\FC_Voxels_AAL_ROI\ROI-corr\LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',this_area,'.mat'));
% 
%         for iRun=1:4
% 
%             iiRun = iiRun + 1;
% 
%             all_runs(iiRun).corr = FC_Voxels.run(iRun).rho_RestingState;
%             all_runs(iiRun).pval = FC_Voxels.run(iRun).pval_RestingState;
% 
%         end
% 
%         clear FC_Voxels
% 
%     end
% 
%     save(strcat('All-Subjects-FC-Voxel-AAL-RestingState-',this_area,'.mat'),'all_runs','-v7.3');

%     clear all_runs

    %% PASSIVE VIEWING

    iiRun = 0;
    for iSubject=1:nSubjects

        disp(int2str(iSubject));
        
        subject = SUBJECT_FOLDER{iSubject};
        
        load(strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\',subject,'\output\FC_Voxels_AAL_ROI\ROI-corr\LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-Passive-',this_area,'.mat'));

        for iRun=1:4

            iiRun = iiRun + 1;

            all_runs(iiRun).corr = FC_Voxels.run(iRun).rho_Passive;
            all_runs(iiRun).pval = FC_Voxels.run(iRun).pval_Passive;

        end

        clear FC_Voxels

    end

    save(strcat('All-Subjects-FC-Voxel-AAL-Passive-',this_area,'.mat'),'all_runs','-v7.3');
    
    clear all_runs
    
    %% TASK
    
    iiRun = 0;
    for iSubject=1:nSubjects

        disp(int2str(iSubject));
        
        subject = SUBJECT_FOLDER{iSubject};
        
        load(strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\',subject,'\output\FC_Voxels_AAL_ROI\ROI-corr\LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-Track-',this_area,'.mat'));

        for iRun=1:4

            iiRun = iiRun + 1;

            all_runs(iiRun).corr = FC_Voxels.run(iRun).rho_Track;
            all_runs(iiRun).pval = FC_Voxels.run(iRun).pval_Track;

        end

        clear FC_Voxels

    end

    save(strcat('All-Subjects-FC-Voxel-AAL-Track-',this_area,'.mat'),'all_runs','-v7.3');
    
    clear all_runs

end