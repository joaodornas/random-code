%% SCRIPT: Run_TVAmodelfitting.m
%
%  Fit TVA to a set of given files (.DAT) in the current directory  
%
%

%% User section
% -----------------

% DATAFILE(S) (.DAT)
    % (NOTE: To run the script from this page, the variable 'Datfile' should 
    % be given. However, the program can also be run from the workspace, given 
    % that a variable 'DAT' is provided. As soon as this variable is present 
    % in the workspace, it will be prioritized over 'Datfile' and the latter 
    % variable will be ignored.) 
    Pathway = '/Volumes/dropbox/_DATA/TVA/Pilot'; 
%    Datfile = {...
             %strcat(Pathway,'/CombiTVA_2346Targets-110-1.dat')};%, ...
             %strcat(Pathway,'/CombiTVA_2346Targets-111-1.dat')};%, ...
             %strcat(Pathway,'/CombiTVA_2346Targets-111-2.dat')};
    %Datfile = {...
             %strcat(Pathway,'/CombiTVA_2346Targets-301-1.dat'), ...
             %strcat(Pathway,'/CombiTVA_2346Targets-302-1.dat'), ...
             %strcat(Pathway,'/CombiTVA_2346Targets-303-1.dat'), ...
             %strcat(Pathway,'/CombiTVA_2346Targets-304-1.dat'), ...
             %strcat(Pathway,'/CombiTVA_2346Targets-305-1.dat')};
%     Datfile = {...
%             strcat(Pathway,'/CombiTVA_2346Targets-302-1.dat'), ...
%             strcat(Pathway,'/CombiTVA_2346Targets-303-1.dat')};
%      Datfile = {...
%             strcat(Pathway,'/CombiTVA_2346Targets-304-1.dat'), ...
%             strcat(Pathway,'/CombiTVA_2346Targets-305-1.dat')};
%               strcat(Pathway,'/CombiTVA_2346Targets-306-1.dat')};

%       Datfile = {strcat(filepath,'/',filename,'.dat')};

     %Datfile = {strcat(Pathway,'/CombiTVA_2346Targets-400-1.dat')};
     %Datfile = {strcat(Pathway,'/CombiTVA_2346Targets-401-1.dat')};
     %Datfile = {strcat(Pathway,'/CombiTVA_2346Targets-402-1.dat')};
     Datfile = {strcat(Pathway,'/CombiTVA_2346Targets-403-1.dat')};

% PROVIDE
    % FOR NEWtvaloader
        % Datafile type (default='STD')
            Type = 'STD';        
        % Unmasked conditions (i.e. condition number) (default=[])
            Unmasked = [];
            
    % FOR NEWtvafit D,L,M,C,P,S,A
        % L (default=1)
            L = 1;        
        % M (default='FREE')
            M = 'FREE';
        % C (default='EXP')
            C = 'EXP';
        % P (default=0)
            P = [];
        % S (default=1)
            S = [];
        % A (default=1)
            A = [];

% SAVE AS
        % Filename
            Filename = 'TVAFITS_Multistability_Pilot.txt';
            
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%% Run
% --------  
Count = 0;
for i=1:length(Datfile)
    if exist(Datfile{i}, 'file')
        Count = Count+1;
    end
end

if Count~=length(Datfile)
    error('One or more datafiles cannot be found');
else
    for i=1:length(Datfile)
        [tvadata{i}, ~] = tvaloader(Datfile{i},Type,Unmasked);
    end

    %matlabpool open
    for i=1:length(Datfile)
        [theta{i},tvamodel{i},tvadata{i}] = tvafit(tvadata{i},L,M,C,P,S,A);
    end
    %matlabpool close

    for i=1:length(Datfile)
        if(i==1) 
            tvalpr(Filename,'',tvadata{i},tvamodel{i},theta{i}); 
        end
        tvalpr(Filename,Datfile{i},tvadata{i},tvamodel{i},theta{i});
        
        disp(strcat('%%% TVA REPORT:',int2str(i)));
        
        tvareport(tvadata{i},tvamodel{i},theta{i});
        
        [t,o,p,c]=tvaplot(tvadata{i},tvamodel{i},theta{i});

        figure; plot(t(1:6),o(1:6),'o'); hold on; plot(t(1:6),p(1:6),'-'); title(strcat('TVA REPORT',int2str(i)));
        
    end
end


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fprintf('\n Run_TVAmodelfitting.m DONE. \n');