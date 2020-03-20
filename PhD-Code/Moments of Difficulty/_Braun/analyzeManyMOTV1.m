
start_ = input('First MOT:');
end_ = input('Last MOT:');

saveAnalyze = true;

prefixdata = 'MOT-10-m-';
prefixanalize = 'MOT-Analyze-10-m-';

for i=start_:end_
    
    disp(strcat('nMOT:',int2str(i)));

    datafile = strcat(prefixdata,int2str(i),'.mat');
    saveAnalyzeDataFile = strcat(prefixanalize,int2str(i),'.mat');
    
    AnalyzeTrajectoriesV1(saveAnalyze,datafile,saveAnalyzeDataFile);
    
end