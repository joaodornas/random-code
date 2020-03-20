%% plot Sparseness distribution

files = dir('Sparseness\*.mat');

nFiles = length(files);

nRes = 14;

nCond = 2;
   
iSpar = 0;
   
for iFile=1:nFiles

   load(strcat('Sparseness\',files(iFile).name));

    for iCond=1:nCond

        iSpar = iSpar + 1;
        sparseness(iSpar) = sparseData(iCond).gallant2002Sparseness;

    end

end
    

f = figure;

histfit(sparseness,10,'kernel');

print(f,'-depsc',strcat('sparseness','.eps'));
    



