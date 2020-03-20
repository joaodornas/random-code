function gallant2005

tic

disp('Load SpikeTrains...');
registros = loadSpikeTrains;

disp('Load Frames...');
stimuli = loadFrames(registros,18);

u = 30;
T = 300;

disp('Load SU...');
stimuli = loadSU(registros,stimuli,u);

disp('Load S...');
stimuli = loadS(registros,stimuli,u);

disp('Load PseudoInverse...');
stimuli = loadPseudoInverse(registros,stimuli);

disp('Load DataSets...');
registros = loadDataSets(registros);

disp('Load Kernel...');
h = loadKernel(registros,stimuli,u);

disp('Prediction...');
results = prediction(registros,h,stimuli,T,u);

disp('Plot Data...');
plotData(results,h);

function registros = loadSpikeTrains
        
    start_time = 500;
    end_time = 9500;
    bin_size = 30;

    registros.labels = { '_nsp006a01_1b-v3' '_nsp006a02_1b-v4' '_nsp006b01_1b-v1' '_nsp007a01_1b-v4' '_nsp007a01_2b-v4' '_nsp007a02_1b-v1' '_nsp007a02_2b-v1' '_nsp007a03_2a-v3' '_nsp008a01_1b-v1' '_nsp008a01_2b-v1' '_nsp008a02_1a-v3' '_nsp008a02_2b-v3' '_nsp008a03_1b-v4' '_nsp008a03_2a-v4' '_nsp009a01_2b-v2' '_nsp009b01_1a-v4' '_nsp010a1_1b-v2' '_nsp010a02_1b-v4' '_nsp011a01_1a-v3' '_nsp011a02_1b-v4' '_nsp012a01_1b-v3' '_nsp012a02_1a-v4' '_nps013a01_1b-v1' '_nsp033a09_3b-v331' 'nsp033a09_1b-v311' 'nsp033a09_1c-v321'};

    %MUA = {'_nsp005a01_1b-v3' '_nsp005a01_2b-v3'};

    % registros(1).label = {'_nsp005a01_1b-v3'};
    % registros(2).label = {'_nsp005a01_2b-v3'};
    
    registros.data(1).label = {'_nsp006a01_1b-v3'};
%     registros.data(2).label = {'_nsp006a02_1b-v4'};
%     registros.data(3).label = {'_nsp006b01_1b-v1'};
%     registros.data(4).label = {'_nsp007a01_1b-v4'};
%     registros.data(5).label = {'_nsp007a02_1b-v1'};
%     registros.data(6).label = {'_nsp007a01_2b-v4'};
%     registros.data(7).label = {'_nsp007a02_2b-v1'};
%     registros.data(8).label = {'_nsp007a03_2a-v3'};
%     registros.data(9).label = {'_nsp008a01_1b-v1'};
%     registros.data(10).label = {'_nsp008a02_1a-v3'};
%     registros.data(11).label = {'_nsp008a03_1b-v4'};
%     registros.data(12).label = {'_nsp008a01_2b-v1'};
%     registros.data(13).label = {'_nsp008a02_2b-v3'};
%     registros.data(14).label = {'_nsp008a03_2a-v4'};
%     registros.data(15).label = {'_nsp009a01_2b-v2'};
%     registros.data(16).label = {'_nsp009b01_1a-v4'}; 
%     registros.data(17).label = {'_nsp010a1_1b-v2'};
%     registros.data(18).label = {'_nsp010a02_1b-v4'};
%     registros.data(19).label = {'_nsp011a01_1a-v3'};
%     registros.data(20).label = {'_nsp011a02_1b-v4'};
%     registros.data(21).label = {'_nsp012a01_1b-v3'};
%     registros.data(22).label = {'_nsp012a02_1a-v4'};
%     registros.data(23).label = {'_nps013a01_1b-v1'};
%     registros.data(24).label = {'nsp033a09_1b-v311'};
%     registros.data(25).label = {'nsp033a09_1c-v321'};
%     registros.data(26).label = {'_nsp033a09_3b-v331'};  

    registros.data(1).videos = 3;
%     registros.data(2).videos = 4;
%     registros.data(3).videos = 1;
%     registros.data(4).videos = 4;
%     registros.data(5).videos = 1;
%     registros.data(6).videos = 4;
%     registros.data(7).videos = 1;
%     registros.data(8).videos = 3;
%     registros.data(9).videos = 1;
%     registros.data(10).videos = 3;
%     registros.data(11).videos = 4;
%     registros.data(12).videos = 1;
%     registros.data(13).videos = 3; 
%     registros.data(14).videos = 4;
%     registros.data(15).videos = 2;
%     registros.data(16).videos = 4;
%     registros.data(17).videos = 2;
%     registros.data(18).videos = 4;
%     registros.data(19).videos = 3;
%     registros.data(20).videos = 4;
%     registros.data(21).videos = 3;
%     registros.data(22).videos = 4;
%     registros.data(23).videos = 1;
%     registros.data(24).videos = 3;
%     registros.data(25).videos = 3;
%     registros.data(26).videos = 3;
    
    registros.data(1).CRF = [604 402];
%     registros.data(2).CRF = [604 402];
%     registros.data(3).CRF = [587 403];
%     registros.data(4).CRF = [679 496];
%     registros.data(5).CRF = [679 496];
%     registros.data(6).CRF = [598 435];
%     registros.data(7).CRF = [598 435];
%     registros.data(8).CRF = [598 435];
%     registros.data(9).CRF = [602 362];
%     registros.data(10).CRF = [602 362];
%     registros.data(11).CRF = [602 362];
%     registros.data(12).CRF = [569 356];
%     registros.data(13).CRF = [569 356];
%     registros.data(14).CRF = [569 356];
%     registros.data(15).CRF = [164 378];
%     registros.data(16).CRF = [587 370];
%     registros.data(17).CRF = [724 443];
%     registros.data(18).CRF = [724 443];
%     registros.data(19).CRF = [689 442];
%     registros.data(20).CRF = [689 442];
%     registros.data(21).CRF = [683 458];
%     registros.data(22).CRF = [683 458];
%     registros.data(23).CRF = [560 408];
%     registros.data(24).CRF = [634 271];
%     registros.data(25).CRF = [634 271];
%     registros.data(26).CRF = [631 273];
    

    for c=1:length(registros.data)
    
        t = 1;
            
            disp(strcat('Registro:',char(registros.data(c).label)));

            %%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Load Spass...');

                Spass = load(strcat(char(registros.data(c).label),'.mat'));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% LOAD CONDITIONS TRIALS label   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Load Conditions Trials label...');

                nConditions = 1;
                
                for i=1:nConditions

                    trials_label(i).label = find(Spass.stimIds == i); 

                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

            %%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Set Resolution...');

                spike_times = Spass.spike_times;

                spike_times = spike_times./ 32000;
                
                clear Spass;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Calculate Spike Rate...');
            
               trials_spikes = spike_times(trials_label(i).label,:);

               nTrials = size(trials_spikes,1);

               for j=1:nTrials

                   trial = trials_spikes(j,:);

                   trial = trial(trial>0);

                   trial = trial(trial>(start_time/1000) & trial<(end_time/1000));

                   nBins = (end_time - start_time)/bin_size;

                   for k=1:nBins

                        spikes = length(trial(trial>=((k-1)*bin_size/1000 + start_time/1000) & trial<(k*bin_size/1000 + start_time/1000)));

                        registros.data(c).Trials(t).rate(k) = spikes / (bin_size/1000);

                   end

                   t = t + 1;

               end
               
               clear trials_spikes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    end
    
    save('registros','registros');

end

function stimuli = loadFrames(registros,pixels)

    width = 1024;
    height = 720;
    nFrames = 300;
    CRFlength = 68;
    
    w = hann(pixels);
    
    wH = zeros(pixels,pixels);
    
    for i=1:pixels
        
        for j=1:pixels
            
            wH(j,i) = w(j);
            
        end
        
    end
    
    w2D = wH.*wH;
    
    clear wH;
    
    for c=1:length(registros.data)
              
      v = registros.data(c).videos;

      disp('Load video...');
      videoin = VideoReader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);

      for f=1:nFrames

          frame = read(videoin,f);

          x = (registros.data(c).CRF(1) - CRFlength):(registros.data(c).CRF(1) + CRFlength);
          y = (registros.data(c).CRF(2) - CRFlength):(registros.data(c).CRF(2) + CRFlength);

          framecrop = frame(x,y);

          frameresize = imresize(framecrop,[pixels pixels]);

          framehann = frameresize.*w2D;

          frameFFT = fft2(framehann);

          frameamp = abs(frameFFT);

          framepower = frameamp^2;

          stimuli.registros.data(c).s(:,:,f) = framepower(:,:);

      end
      
      clear videoin;
        
    end
    
    save('stimuli','stimuli');

end

function stimuli = loadSU(registros,stimuli,u)

    for c=1:length(registros.data)
        
        T = size(stimuli.registros.data(c).s,3);
        
        disp('Calculate SU...');
        for U=0:(u-1)  
            
            for i=1:T;

                n = 1;
                for j=1:size(stimuli.registros.data(c).s,1)

                    for k=1:size(stimuli.registros.data(c).s,2)
                        
                        lag = i - U;
                        
                        if lag > 0
                            
                            stimuli.registros.data(c).SU(i,n,U+1) = stimuli.registros.data(c).s(j,k,lag);
                            
                        else 
                            
                            stimuli.registros.data(c).SU(i,n,U+1) = 0;
                            
                        end
                        
                        n = n + 1;
                        
                    end

                end

            end
            
        end
        
        disp(strcat('Size SU:', num2str(size(stimuli.registros.data(1).SU))));
        
        clear SU;
        
    end
    
    save('stimuli','stimuli');
    
end
    
function stimuli = loadS(registros,stimuli,u)
        
    for c=1:length(registros.data)
        
        N = size(stimuli.registros.data(c).s,1)*size(stimuli.registros.data(c).s,2);

        T = size(stimuli.registros.data(c).s,3);

        disp('Calculate S...');

        x = 1;
        S = zeros(T*u,1:N);
        for U=0:(u-1)
            
            for t=1:T
                
                line = zeros(N,1);

                line = stimuli.registros.data(c).SU(t,:,(U+1)); 

                S(x,1:N) = line;

                x = x + 1;
                
            end
            
        end

        disp(strcat('Size S:',num2str(size(S))));
        
        Slinha = transpose(S);
        
        stimuli.registros.data(c).S = Slinha;

        clear S;
        
    end
    
    save('stimuli','stimuli');

end

function stimuli = loadPseudoInverse(registros,stimuli)
        
    for c=1:length(registros.data)
        
        T = size(stimuli.registros.data(c).s,3);
        
        beta = linspace(1/100000,1/10,30);

        stimuli.registros.data(c).beta = beta;

        disp('Calculate Pseudo-Inverse...');
        for i=1:length(beta)

           pseudo(:,:,i) = pinv(((stimuli.registros.data(c).S.*stimuli.registros.data(c).S)./T),beta(i));      

        end

        disp(strcat('Size Pseudo-Inverse:',num2str(size(pseudo))));
        
        stimuli.registros.data(c).pseudo = pseudo;
    
    end
    
    save('stimuli','stimuli');
    
end

function registros = loadDataSets(registros)
        
    for c=1:length(registros.data)
        
        nTrials = length(registros.data(c).Trials);

        nDataSets = 20;

        jackknife = round((5/100)*nTrials);
        
        for i=1:nDataSets
            
            for j=1:nTrials
                
                comb(i,j) = randi(nTrials);
            
            end
            
        end

        disp('Calculate DataSets...');
        for i=1:nDataSets

            for j=1:nTrials
                   
                dataSets(i).trial(j).rate = registros.data(c).Trials(comb(i,j)).rate; 

            end

            dataSets(i).trial = dataSets(i).trial(1:end-jackknife);

            dataSets(i).validation = dataSets(i).trial(end-jackknife:end);

        end
        
        registros.data(c).dataSets = dataSets;
        
        clear dataSets;
    
    end
    
    save('registros','registros');
    
end

function h = loadKernel(registros,stimuli,u)
        
    for c=1:length(registros.data)
        
        beta = length(stimuli.registros.data(c).beta);
        
        nDataSets = length(registros.data(c).dataSets);
        
        T = size(stimuli.registros.data(c).s,3);
        
        disp('Calculate Kernel...');
        for i=1:length(beta)

            for j=1:nDataSets

                nTrials = length(registros.data(c).dataSets(j).trial);

                rp = zeros(T,1);

                for k=1:nTrials

                    for l=1:T

                        rp(l) = rp(l) + registros.data(c).dataSets(j).trial(k).rate(l);

                    end

                end

                nValidation = length(registros.data(c).dataSets(j).validation);

                rv = zeros(T,1);

                for k=1:nValidation

                    for l=1:T

                        rv(l) = rv(l) + registros.data(c).dataSets(j).validation(k).rate(l);

                    end

                end

                rp = rp./nTrials;

                rv = rv./nValidation;

                h.registros.data(c).beta(i).dataSet(j).rp = rp;

                h.registros.data(c).beta(i).dataSets(j).rv = rv;

                h.registros.data(c).beta(i).dataSet(j).dataSets = registros.data(c).dataSets(j);

                h.registros.data(c).beta(i).dataSet(j).betaValue = beta(i);
                
                disp(strcat('Size rv:', num2str(size(rv))));
                disp(strcat('Size rp:', num2str(size(rp))));
                
                cs = stimuli.registros.data(c).pseudo(:,:,i)*stimuli.registros.data(c).S;
                
                disp(strcat('Size cs:',num2str(size(cs))));
                
                rph = [];
                for r=1:u
                    
                    rph = [rph; rp];
                    
                end
                
                disp(strcat('Size rph:', num2str(size(rph))));
                
                h.registros.data(c).beta(i).dataSet(j).kernel = cs*rph;

            end

            meanKernel = 0;
            for j=1:nDataSets

                meanKernel = meanKernel + h.registros.data(c).beta(i).dataSet(j).kernel;

            end

            meanKernel = meanKernel./nDataSets;

            stdKernel = 0;
            for j=1:nDataSets

                stdKernel = stdKernel + (h.registros.data(c).beta(i).dataSet(j).kernel - meanKernel).^2;

            end

            stdKernel = sqrt(stdKernel./nDataSets);

            h.registros.data(c).beta(i).meanKernel = meanKernel;

            h.registros.data(c).beta(i).stdKernel = stdKernel;

            gammas = linspace((8/10),2,7);

            H = [];
            for g=1:length(gammas)

                rectification = 1 - gammas(g).*stdKernel.^2/meanKernel.^2;

                if rectification < 0

                    rectification = 0;

                end
                
                H(g).kernel = (meanKernel.')*sqrt(rectification);

                H(g).gamma = gammas(g);

            end

            h.registros.data(c).beta(i).H = H;

            h.registros.data(c).gammas = gammas;

        end
   
    end
    
    save('h','h');
    
end

function results = prediction(registros,h,stimuli,T,u)

    for c=1:length(registros.data)
        
       nBeta = length(h.registros.data(c).beta);
       
       disp('Extract Data...');
       
       disp(strcat('Size S:',num2str(size(stimuli.registros.data(c).S))));
       disp(strcat('Size kernel:',num2str(size(h.registros.data(c).beta(1).H(1).kernel))));
               
       for i=1:nBeta
          
           nGamma = length(h.registros.data(c).gammas);
           
           for g=1:nGamma
               
               n = size(stimuli.registros.data(c).S.',2);
               
               kernel = [];
               for z=1:n
                   
                   kernel = [kernel; h.registros.data(c).beta(i).H(g).kernel];
                   
               end
               
               SH = stimuli.registros.data(c).S.'*kernel;
               
               disp(strcat('Size SH:',num2str(size(SH))));
               
               teta = linspace(min(SH),max(SH),1);
               
%                residual = zeros(length(SH));
               
%                for k=1:length(SH)
%                   
%                    residual(k) = SH(k) - mean(SH);
%                    
%                end
               
               for t=1:length(teta)
                  
                   %response = SH - teta(t) + residual;
                   
                   responseh = SH - teta(t);
                   
                   disp(strcat('Size responseh:',num2str(size(responseh))));
                   
                   for l=1:length(responseh)
                      
                       if responseh(l) < 0
                           
                           responseh(l) = 0;
                           
                       end
                       
                   end
                   
                   response = zeros(1,T);
                   
                   disp(strcat('Size response:',num2str(size(response))));
                   
                   for t=1:T*u
                       
                       for b=0:(u-1)
                           
                           %disp(strcat('Size responseh(t,((n*u)+1):((n*u)+u)):',num2str(size(responseh(t,((n*u)+1):((n*u)+u))))));
                           
                           response = response + responseh(t,((b*T)+1):((b*T)+T));
                       
                       end
                       
                   end
                   
                   response = response ./ (u*T*u);
                                      
                   for y=1:length(h.registros.data(c).beta(i).dataSets)
                               
                       disp(strcat('Size rv:', num2str(size(h.registros.data(c).beta(i).dataSets(y).rv))));
                       
                       [rho,pval] = corr(transpose(response),h.registros.data(c).beta(i).dataSets(y).rv,'type','pearson');
                       
                       results.registros.data(c).beta(i).gamma(g).teta(t).tetaValue = teta(t);
                       
                       results.registros.data(c).beta(i).gamma(g).teta(t).rho = rho;
                       
                       results.registros.data(c).beta(i).gamma(g).teta(t).pval = pval;
                       
                       results.registros.data(c).beta(i).gamma(g).response = response;
                                              
                   end
                   
               end
               
           end
           
       end        
        
    end    
        
    save('results','results');
    
end

function plotData(results,h)

    for c=1:length(results.registros.data)
        
        rho = 0;
        pval = 0;
        
        rpre = [];
        rv = [];
        
        for i=1:length(results.registros.data(c).beta)
            
            for g=1:length(results.registros.data(c).beta(i).gamma)
                
                for t=1:length(results.registros.data(c).beta(i).gamma(g).teta)
                    
                    if rho < results.registros.data(c).beta(i).gamma(g).teta(t).rho
                        
                        rho = results.registros.data(c).beta(i).gamma(g).teta(t).rho;
                        pval = results.registros.data(c).beta(i).gamma(g).teta(t).pval;
                        
                        rpre = results.registros.data(c).beta(i).gamma(g).response;
                        
                        for j=1:length(h.registros.data(c).beta(i).dataSets)
                            
                            rv = rv + h.registros.data(c).beta(i).dataSets(j).rv;
                            
                        end
                        
                        rv = rv./length(h.registros.data(c).beta(i).dataSets);
                
                    end
                
                end
            
            end
    
        end
        
        f = figure;
        
        plot(1:length(rpre),rpre,'r');
        hold on;
        
        plot(1:length(rv),rv,'b');
        
        print(f,'-depsc',strcat('registro-',int2str(c)));
        
    end

end

toc

end

