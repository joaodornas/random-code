%% plot Reliability distribution

% % % files = dir('Reliability\*.mat');
% % % 
% % % nFiles = length(files);
% % % 
% % % nRes = 14;
% % % 
% % % for iRes=1:nRes
% % %    
% % %    %nCond = length(reliabilityPerCondition);
% % %    nCond = 2;
% % %    
% % %    iReal = 0;
% % %    
% % %    for iFile=1:nFiles
% % %        
% % %        load(strcat('Reliability\',files(iFile).name));
% % %        
% % %         for iCond=1:nCond
% % %        
% % %             iReal = iReal + 1;
% % %             real(iRes,iReal) = reliabilityPerCondition(iCond).resolution(iRes).reliability;
% % %        
% % %         end
% % %         
% % %    end
% % %     
% % %    clear reliabilityPerCondition
% % %     
% % % end
% % % 
% % % for iRes=1:nRes
% % %    
% % %     f = figure;
% % %     
% % %     histfit(real(iRes,:),10,'kernel');
% % %     
% % %     xlim([0 1]);
% % %     
% % %     print(f,'-depsc',strcat('resolution-',int2str(iRes),'.eps'));
% % %     
% % %     close all
% % %     
% % % end
% % % 
% % % clear real

%% Por Protocol and Condition

% % % protocol = {'VideosFB','VideosCC'};
% % % nProtocol = length(protocol);
% % % 
% % % for iProtocol=1:nProtocol
% % %     
% % %     [pathToData,dataFolders] = loadPathToData(protocol{iProtocol});
% % % 
% % %     nData = length(dataFolders);
% % % 
% % %     for iData=1:nData
% % % 
% % %         dataFile = dataFolders{iData};
% % % 
% % %         idx = strfind(dataFile,'\');
% % % 
% % %         dataFile = dataFile(idx(end)+1:end);
% % % 
% % %         load(strcat('Reliability\',dataFile,'-reliability.mat'));
% % % 
% % %         for iCond=1:nCond
% % % 
% % %             for iRes=1:nRes
% % % 
% % %                 real(iCond,iRes,iData) = reliabilityPerCondition(iCond).resolution(iRes).reliability;
% % % 
% % %             end
% % % 
% % %         end
% % % 
% % %     end
% % % 
% % %     for iCond=1:nCond
% % % 
% % %         for iRes=1:nRes
% % % 
% % %             f = figure;
% % % 
% % %             histfit(squeeze(real(iCond,iRes,:)),10,'kernel');
% % % 
% % %             xlim([0 1]);
% % % 
% % %             print(f,'-depsc',strcat(protocol{iProtocol},'-resolution-',int2str(iRes),'-iCond-',int2str(iCond),'.eps'));
% % % 
% % %             close all
% % % 
% % %         end
% % % 
% % %     end
% % % 
% % % end

%% Por Protocol

protocol = {'VideosFB','VideosCC'};
nProtocol = length(protocol);

nRes = 14;
nCond = 2;

for iRes=1:nRes
    
    disp(int2str(iRes));

        [pathToData,dataFolders] = loadPathToData('VideosCC');

        nData = length(dataFolders);
        
        iD = 0;

        for iData=1:nData

            dataFile = dataFolders{iData};

            idx = strfind(dataFile,'\');

            dataFile = dataFile(idx(end)+1:end);

            load(strcat('Reliability\',dataFile,'-reliability.mat'));

            for iCond=1:nCond
                
                iD = iD + 1;

                real(iRes,iD) = reliabilityPerCondition(iCond).resolution(iRes).reliability;

            end

        end
        
end

for i=1:nRes res(i) = reliabilityPerCondition(1).resolution(i).resolution; end

errorbar(1:nRes,nanmean(real,2),nanstd(real,0,2));

set(gca,'xtick',1:nRes);
set(gca,'xticklabel',res);
set(gca,'xlim',[1 nRes]);

