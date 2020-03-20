
velocities = [20 25 30 35 40 45 50];
nMOT = 10;

saveAnalyze = true;

prefixdata = 'MOT-';
prefixanalize = 'MOT-Analyze-V2-11-m-Collinear-';

for iVelocity=velocities
    
    disp(strcat('nMOT:',int2str(iVelocity)));

    for iMOT=1:nMOT
        
        datafile = strcat(prefixdata,int2str(iVelocity),'-',int2str(iMOT),'.mat');
        saveAnalyzeDataFile = strcat(prefixanalize,int2str(iVelocity),'-',int2str(iMOT),'.mat');
    
        AnalyzeCollinearityV14to2(saveAnalyze,datafile,saveAnalyzeDataFile);
        
    end
    
end