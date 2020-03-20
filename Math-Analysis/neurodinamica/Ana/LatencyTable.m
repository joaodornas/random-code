cells = {'gas004a01_3a' 'gas004a02_1c' 'gas005a01_1b' 'gas005c01_1b' 'gas006a02_1b' 'gas008a01_1c' 'gas009a02_2b' 'gas009a02_2c' 'gas009a04_3a' 'gas010a02_2b' 'gas011b02_1b' 'gas012a02_1b' 'gas013a02_2c' 'gas016a02_3a' 'gas020b02_1b' 'gas024a03_2B' 'gas024a03_2C' 'gas025a02_1B' 'gas025a02_1C' 'gas026a03_3B' 'gas031a02_1B' 'gas034a02_1B' 'gas036a03_2B' 'gas036a03_2C' 'gas037a03_1B' 'gas038b02_1B' 'gas039a02_1C' 'gas040a01_2B' 'gas043b02_2B' 'gas046b02_1B' 'gas049a02_2C' 'gas051a02_1B' 'gas051d02_1B' 'gas053b02_1A' 'gas053b03_2B' 'gas054b01_2B' 'gas055a02_1B' 'gas056a02_1B' 'gas057a02_1B' 'gas057a02_1C' 'gas057b03_2B' 'gas058b02_1A' 'gas059a02_1B' 'gas061a03_2B' 'gas061b02_2B' 'gas061b03_1A' 'gas061b03_1C' 'gas062b02_1B' 'gas063a02_1B' 'gas066b01_1B' 'gas068b02_1A' 'gas071a02_2B' 'gas071a03_2A' 'gas072a05_1B' 'gas075a02_3A' 'gas075b02_2B' 'gas076a02_1B' 'gas081b02_1B' 'gas083a02_1B' 'gas083a04_3B' 'gas084a04_3b' 'gas087a02_3B' 'gas091a02_3a' 'gas091b02_3b' 'gas093a02_2B' 'gas094a02_1B' 'gas095a02_1B' 'gas096a02_1B' 'gas109a03_2b' 'gas110a03_2a'};
startCondition = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 16 16 16 16 18 18 18 18 1 18 1 1 18 18 18 18 1 1];
endCondition = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 30 30 30 30 34 34 34 34 10 34 17 17 34 34 34 34 16 16];
totalConditions = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 30 30 30 30 34 34 34 34 10 34 17 17 34 34 34 34 32 32];

ncells=numel(cells);
for k=1:ncells
    
    disp(['computing cell #' num2str(k) '...'])
    
    cellname=char(cells(k));
    
    filename=['/Users/altmaia/Flocke/Doutorado/ANALISE/PSTHKernel_withLatency_wholeTime/' cellname '.mat']; 
    
    load(filename);
    
    LatencyTable{k,1}=cellname;
    
    %nconditions=numel(startCondition(k):endCondition(k));
    
    s = 1;
    for i=startCondition(k):endCondition(k)
       
        histData.condition(i).kernel.latencyPoints = histData.condition(i).kernel.latencyPoints - 1000/1000;
       
        LatencyTable{k,s+1}=histData.condition(i).kernel.latencyPoints;
        
        s = s + 1;
        
    end
    
    clear histData;
    
end