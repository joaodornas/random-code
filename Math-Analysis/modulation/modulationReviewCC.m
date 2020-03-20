
function modulationReviewCC(allCells)

tic

disp('BEGIN');

excluidos = { '_nsp015a01_2a-v1' '_nsp015a02_2b-v3' '_nsp016a01_2b-v2' '_nsp016a02_2b-v4' '_nsp017a01_1a-v3' '_nsp017b01_2b-v2' '_nsp018b01_1b-v3' };
registros = { '_nsp020a01_1b-v5' '_nsp020a02_1b-v6' '_nsp021a01_1a-v5' '_nsp021a02_1a-v6' '_nsp022a01_1b-v5' '_nsp022a02_1b-v6' '_nsp023a01_1a-v1' '_nsp023a02_1a-v4' '_nsp024a01_1a-v6' '_nsp024a02_1a-v3' '_nsp025a01_1b-v2' '_nsp025a02_1a-v1' '_nsp026a01_3b-v8' '_nsp026a02_3a-v11' '_nsp027a03_1a-v10' '_nsp027a04_1b-v9' '_nsp028a03_1a-v6' '_nsp028a04_1b-v8' '_nsp029a03_1a-v8' '_nsp029a04_1b-v3' '_nsp030a03_2a-v5' '_nsp030a04_2b-v8' '_nsp031b03_1b-v3' '_nsp031b04_1b-v4' '_nsp032a03_1a-v8' '_nsp032a04_1a-v2' '_nsp033a04_1b-v821' '_nsp033a04_1c-v811' '_nsp033a04_3a-v831' '_nsp033a05_1b-v511' '_nsp033a05_1c-v521' '_nsp033a05_3a-v531' '_nsp033a06_1b-v812' '_nsp033a06_1c-v822' '_nsp033a07_1b-v512' '_nsp033a07_3b-v532' '_nsp034a04_2b-v8' '_nsp034a05_2b-v11' '_nsp035a03_2b-v8' '_nsp035a04_2b-v1'};

start_time = 500;
end_time = 9500;
bin_size = 30;


%%%%%%%%%%%%%%%%%%%%%%% PROCESSA TODOS OS REGISTROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if allCells == 1

    for c=1:length(registros)
    
        %%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Load Spass...');

            Spass = load(strcat(char(registros(c)),'.mat'));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Load Conditions Trials Labels...');

        nConditions = max(Spass.stimIds);
        
            for i=1:nConditions

                trials_label(i).label = find(Spass.stimIds == i); 

            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        %%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Set Resolution...');

            spike_times = Spass.spike_times;

            spike_times = spike_times./ 32000;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Begin Conditions Loop...');

            for i=1:nConditions

               disp(strcat('Begin Condition . ',int2str(i))); 

               trials_spikes = spike_times(trials_label(i).label,:);
               
               nTrials = size(trials_spikes,1);

               for j=1:nTrials

                   trial = trials_spikes(j,:);

                   trial = trial(trial>0);

                   trial = trial(trial>(start_time/1000) & trial<(end_time/1000));

                   nBins = (end_time - start_time)/bin_size;

                   for k=1:nBins

                        spikes = length(trial(trial>=((k-1)*bin_size/1000 + start_time/1000) & trial<(k*bin_size/1000 + start_time/1000)));

                        condition(i).allTrials(j).binRateTrial(k) = spikes/(bin_size/1000);
                        condition(i).allTrials(j).binRateTrialSq(k) = ( spikes/(bin_size/1000) )^2;

                   end

                   condition(i).meanFRTrial(j) = mean(condition(i).allTrials(j).binRateTrial);
                   condition(i).stdFRTrial(j) = std(condition(i).allTrials(j).binRateTrial); 
                   
                   condition(i).sparseness.Trial(j).gallant2002SelectivityIndex = ( 1 - ( ( sum(condition(i).allTrials(j).binRateTrial) )^2 / ( nBins * sum(condition(i).allTrials(j).binRateTrialSq) ) ) ) / ( 1 - ( 1/nBins ) ) ;
                   condition(i).sparseness.Trial(j).gallant2002Sparseness = ( 1 - ( mean(condition(i).allTrials(j).binRateTrial)^2 / ( mean(condition(i).allTrials(j).binRateTrial)^2 + std(condition(i).allTrials(j).binRateTrial)^2 ) ) / ( 1 - ( 1 / nBins ) ) ) ;
                   condition(i).sparseness.Trial(j).kurtosis = kurtosis(condition(i).allTrials(j).binRateTrial) - 3; 
                   
               end

               spikes_vector = reshape(trials_spikes.',[],1);

               spikes_vector = spikes_vector.';

               spikes_vector = sort(spikes_vector);

               spikes_vector = spikes_vector(spikes_vector>0);

               spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));

               nBins = (end_time - start_time)/bin_size;
               
                   for k=1:nBins

                        spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));

                        condition(i).meanAllTrials.meanBinRate(k) = spikes/((nTrials)*(bin_size/1000));
                        condition(i).meanAllTrials.meanBinRateSq(k) = ( spikes/((nTrials)*(bin_size/1000)) )^2;

                   end

                   condition(i).meanFRAllTrials = mean(condition(i).meanAllTrials.meanBinRate);
                   condition(i).stdFRAllTrials = std(condition(i).meanAllTrials.meanBinRate);
                   
                   condition(i).sparseness.gallant2002SelectivityIndex = ( 1 - ( ( sum(condition(i).meanAllTrials.meanBinRate) )^2 / ( nBins * sum(condition(i).meanAllTrials.meanBinRateSq) ) ) ) / ( 1 - ( 1/nBins ) ) ;
                   condition(i).sparseness.gallant2002Sparseness = ( 1 - ( mean(condition(i).meanAllTrials.meanBinRate)^2 / ( mean(condition(i).meanAllTrials.meanBinRate)^2 + std(condition(i).meanAllTrials.meanBinRate)^2 ) ) / ( 1 - ( 1 / nBins ) ) ) ;
                   condition(i).sparseness.kurtosis = kurtosis(condition(i).meanAllTrials.meanBinRate) - 3;               

            end
  
            nCRF = 2;
            CRF = 1;
          
            for k=1:nBins
                
               diff.nCRFxCRF.difBins(k) = condition(nCRF).meanAllTrials.meanBinRate(k) - condition(CRF).meanAllTrials.meanBinRate(k);
               
            end
            
            diff.nCRFxCRF.meanDifBins = mean(diff.nCRFxCRF.difBins);
            diff.nCRFxCRF.stdDifBins = std(diff.nCRFxCRF.difBins); 
               
            sparseGain.nCRFxCRF = ( condition(nCRF).sparseness.gallant2002SelectivityIndex / condition(CRF).sparseness.gallant2002SelectivityIndex )*100 - 100;
               
            if nConditions > 2

                for k=1:nBins
                    
                    diff.ExtraNCRFxCRF.difBins(k) = condition(3).meanAllTrials.meanBinRate(k) - condition(1).meanAllTrials.meanBinRate(k);
                    
                end
                
                diff.ExtraNCRFxCRF.meanDifBins = mean(diff.nCRFxCRF.difBins);
                diff.ExtraNCRFxCRF.stdDifBins = std(diff.nCRFxCRF.difBins); 

                sparseGain.ExtraNCRFxCRF = ( condition(3).sparseness.gallant2002SelectivityIndex / condition(1).sparseness.gallant2002SelectivityIndex )*100 - 100;
               
            end
            
            modulation.sparseGain = sparseGain;
            modulation.conditions = condition;
            modulation.diff = diff;
            
            filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1');

            mkdir(filepath,'ModulationReview');

            filepath = strcat(filepath,'/','ModulationReview','/',char(registros(c)));

            filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-modulationReview');

            save(filepath,'modulation');
            
            clear condition;
            clear sparseGain;
            clear diff;
            clear modulation;
            
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% ANALISA TODOS OS REGISTROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = 1;
sparseSig = 0;
for c=1:length(registros)
    
    folder = '/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/ModulationReview/';
   
    modulation = load(strcat(folder,char(registros(c)),'-bin_size-',int2str(bin_size),'-modulationReview.mat'));
    
    resultadoAll(s) = modulation.modulation.diff.nCRFxCRF.meanDifBins;
    
    sparse(s) = modulation.modulation.sparseGain.nCRFxCRF;
    
    nTrials = min(length(modulation.modulation.conditions(1).sparseness.Trial),length(modulation.modulation.conditions(2).sparseness.Trial));
    
    if length(modulation.modulation.conditions) > 2
        
        nTrials = min(nTrials,length(modulation.modulation.conditions(3).sparseness.Trial));
        
    end
    
    for j=1:nTrials
        
        reg(c).sparseTrialsCRF(j) = modulation.modulation.conditions(1).sparseness.Trial(j).gallant2002SelectivityIndex;
        reg(c).sparseTrialsnCRF(j) = modulation.modulation.conditions(2).sparseness.Trial(j).gallant2002SelectivityIndex;
        
    end
    
    [p,h] = signrank(reg(c).sparseTrialsCRF,reg(c).sparseTrialsnCRF);
    
    reg(c).nCRFxCRFp = p;
    reg(c).nCRFxCRFh = h;
    
    if (h == 1 && p < 0.05), sparseSig = sparseSig + 1; end
        
    s = s + 1;
        
    if length(modulation.modulation.conditions) > 2
    
        resultadoAll(s) = modulation.modulation.diff.ExtraNCRFxCRF.meanDifBins;
        
        sparse(s) = modulation.modulation.sparseGain.ExtraNCRFxCRF;
        
        for j=1:nTrials
        
            reg(c).sparseTrialsExtraNCRF(j) = modulation.modulation.conditions(3).sparseness.Trial(j).gallant2002SelectivityIndex;
        
        end
    
        [p,h] = signrank(reg(c).sparseTrialsCRF,reg(c).sparseTrialsExtraNCRF);
        
        reg(c).ExtraNCRFxCRFp = p;
        reg(c).ExtraNCRFxCRFh = h;
        
        if (h == 1 && p < 0.05) , sparseSig = sparseSig + 1; end
        
        s = s + 1;    
        
    end
    
end

NsupressaoAll = length(resultadoAll(resultadoAll < 0))
NfacilitacaoAll = length(resultadoAll(resultadoAll > 0))

resultados.SupFacAll = resultadoAll;

resultados.nSupressaoAll = NsupressaoAll;
resultados.nFacilitacaoAll = NfacilitacaoAll;

supressaoAll = resultadoAll(resultadoAll < 0);
facilitacaoAll = resultadoAll(resultadoAll > 0);

[p, h] = ranksum(abs(supressaoAll),facilitacaoAll)

resultados.SupFacRanksumP = p;
resultados.SupFacRanksumH = h;

idxSupressao = find(resultadoAll(resultadoAll < 0));
idxFacilitacao = find(resultadoAll(resultadoAll > 0));

sparseSupressao = sparse(idxSupressao);
sparseFacilitacao = sparse(idxFacilitacao);

resultados.sparseness = sparse;
resultados.sparseSupressao = sparseSupressao;
resultados.sparseFacilitacao = sparseFacilitacao;

meanSparseSupressao = mean(sparseSupressao)
meanSparseFacilitacao = mean(sparseFacilitacao)

resultados.meanSparseSupressao = meanSparseSupressao;
resultados.meanSparseFacilitacao = meanSparseFacilitacao;

[p, h] = ranksum(sparseSupressao,sparseFacilitacao)

resultados.sparseSupFacRanksumP = p;
resultados.sparseSupFacRanksumH = h;

[r, p] = pearson(sparseSupressao,supressaoAll);

resultados.corrSparseSupR = r;
resultados.corrSparseSupP = p;

[r, p] = pearson(sparseFacilitacao,facilitacaoAll);

resultados.corrSparseFacR = r;
resultados.corrSparseFacP = p;

resultados.sparseSig = sparseSig;
resultados.reg = reg;

save(strcat(folder,'resultados-modulation-review-CC'),'resultados');

f = figure;

plot(sparseSupressao,supressaoAll,'ro');
g = fittype('a*x + b','coeff',{'a','b'});
[X, gofReverso] = fit(sparseSupressao.',supressaoAll.',g);
hold on;
plot(X);
xlim([0 max(sparseSupressao)]);
ylim([min(supressaoAll) 0]);
xlabel('Sparseness (supressao)');
ylabel('Supressao');
text(max(sparseSupressao)-20,max(supressaoAll),strcat('r =',int2str(round(resultados.corrSparseSupR*100)),'/','p =',num2str(resultados.corrSparseSupP)));

print(f,'-depsc',strcat(folder,'corr-sparseSupressao'));

u = figure;

plot(sparseFacilitacao,facilitacaoAll,'ro');
g = fittype('a*x + b','coeff',{'a','b'});
[X, gofReverso] = fit(sparseFacilitacao.',facilitacaoAll.',g);
hold on;
plot(X);
xlim([0 max(sparseFacilitacao)]);
ylim([0 max(facilitacaoAll)]);
xlabel('Sparseness (facilitacao)');
ylabel('Facilitacao');
text(max(sparseFacilitacao)-20,max(facilitacaoAll),strcat('r =',int2str(round(resultados.corrSparseFacR*100)),'/','p =',num2str(resultados.corrSparseFacP)));

print(u,'-depsc',strcat(folder,'corr-sparseFacilitacao'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc

end
       