function resolutionAlongTrain(stimulus_start_time,intervals)
 
tic

date = {'12-08-20' '12-08-27' '12-08-27' '12-08-29' '12-08-29' '12-08-29' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-09-03' '12-09-03' '12-09-04' '12-09-04' '12-09-05' '12-09-05' '12-09-06' '12-09-06' '12-10-16' '12-12-17' '12-12-17' '12-12-17'};

sitio = {int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(2) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(2) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1)};

channel = {'E1' 'E1' 'E2' 'E1' 'E1' 'E1' 'E1' 'E3' 'E1' 'E3' 'E3' 'E1' 'E3' 'E1' 'E3' 'E1' 'E3' 'E3' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E3'};

video = {int2str(3) int2str(3) int2str(3) int2str(3) int2str(4) int2str(1) int2str(4) int2str(4) int2str(1) int2str(1) int2str(3) int2str(1) int2str(1) int2str(3) int2str(3) int2str(4) int2str(4) int2str(2) int2str(4) int2str(4) int2str(2) int2str(3) int2str(4) int2str(3) int2str(4) int2str(1) int2str(311) int2str(321) int2str(331)};

label = {'nsp003a01_1b' 'nsp005a01_1b' 'nsp005a01_2b' 'nsp006a01_1b' 'nsp006a02_1b' 'nsp006b01_1b' 'nsp007a01_1b' 'nsp007a01_2b' 'nsp007a02_1b' 'nsp007a02_2b' 'nsp007a03_2a' 'nsp008a01_1b' 'nsp008a01_2b' 'nsp008a02_1a' 'nsp008a02_2b' 'nsp008a03_1b' 'nsp008a03_2a' 'nsp009a01_2b' 'nsp009b01_1a' 'nsp010a02_1b' 'nsp010a1_1b' 'nsp011a01_1a' 'nsp011a02_1b' 'nsp012a01_1b' 'nsp012a02_1a' 'nps013a01_1b' 'nsp033a09_1b' 'nsp033a09_1c' 'nsp033a09_3b'};

end_time = [4000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000];

stimulus_end_time = end_time - 500;

    for c=1:length(date)
   
        begin = round( (stimulus_start_time/1000) / (intervals/1000) ) - 1 ;

        finish = round( (stimulus_end_time(c)/1000) / (intervals/1000) ) -1 ;

        origem(c).label = strcat(char(date(c)),'-','sitio',char(sitio(c)),'-',char(channel(c)),'-','v',char(video(c)));
        
        disp(strcat(origem(c).label))
        
        s = 1;
        for i=begin:(finish)
           
            origem(c).Protocol(s).metric = load(strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',char(date(c)),'/sitio',char(sitio(c)),'/',char(channel(c)),'/v',char(video(c)),'/full-movie/trecho-invertido/v',char(video(c)),'-',char(date(c)),'-sitio',char(sitio(c)),'-',char(channel(c)),'-',int2str(intervals*i),'-',int2str(intervals*i+intervals),'-90-bin-size-1-trecho-invertido.mat'));

            s = s + 1;

        end

        s = s - 1;

        Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];
    
        cell(c).label = origem(c).label;
        
        cell(c).nTimingSignificativo = 0;
        cell(c).nRateSignificativo = 0;
        
        g = 1;
        h = 1;
        
            for i=1:length(origem(c).Protocol)

                    cell(c).Qidx(i) = origem(c).Protocol(i).metric.metric_analysis.max_info.max_info_tpmc_idx_T;
                    cell(c).HTiming(i) = origem(c).Protocol(i).metric.metric_analysis.max_info.max_info_tpmc_T;
                    cell(c).HCount(i) = origem(c).Protocol(i).metric.metric_analysis.max_info.Hcount_tpmc_info_T;
                    cell(c).HBias_T(i) = origem(c).Protocol(i).metric.metric_analysis.max_info.HBias_T(cell(c).Qidx(i));
                    cell(c).HBias_std_T(i) = origem(c).Protocol(i).metric.metric_analysis.max_info.HBias_std_T(cell(c).Qidx(i));

                if cell(c).HTiming(i) > cell(c).HCount(i)

                    cell(c).Timing(i) = 1;
                    cell(c).Rate(i) = 0;
                    cell(c).entropy(i) = cell(c).HTiming(i);
                    
                    if cell(c).HTiming(i) > (cell(c).HBias_T(i) + 2*cell(c).HBias_std_T(i))

                        cell(c).significativo(i) = 1;
                        
                        cell(c).nTimingSignificativo = cell(c).nTimingSignificativo + 1;
                        
                        cell(c).ResolutionTimingSignificativo(g) = cell(c).Qidx(i);
                        
                        cell(c).EntropyTimingSignificativo(g) = cell(c).HTiming(i);
                        
                        g = g + 1;

                    else

                        cell(c).significativo(i) = 0;

                    end

                else

                    cell(c).Timing(i) = 0;
                    cell(c).Rate(i) = 1;
                    cell(c).entropy(i) = cell(c).HCount(i);
                    
                    if cell(c).HCount(i) > (cell(c).HBias_T(i) + 2*cell(c).HBias_std_T(i))

                        cell(c).significativo(i) = 1;
                        
                        cell(c).nRateSignificativo = cell(c).nRateSignificativo + 1;
                        
                        cell(c).EntropyRateSignificativo(h) = cell(c).HCount(i);
                        
                        h = h + 1;

                    else

                        cell(c).significativo(i) = 0;

                    end

                end

            end

            cell(c).meanEntropy = mean(cell(c).entropy);
            cell(c).stdEntropy = std(cell(c).entropy);
            cell(c).sumEntropy = sum(cell(c).entropy);
            
            cell(c).meanEntropySignificativo = mean([cell(c).EntropyTimingSignificativo cell(c).EntropyRateSignificativo]);
            cell(c).stdEntropySignificativo = std([cell(c).EntropyTimingSignificativo cell(c).EntropyRateSignificativo]);
            
            cell(c).nTimingSignificativo = cell(c).nTimingSignificativo / s ;
            cell(c).nRateSignificativo = cell(c).nRateSignificativo / s ;
            
            Qidx = cell(c).Qidx./max(cell(c).Qidx);
            entropy = cell(c).entropy./max(cell(c).entropy);
            
            cell(c).mostFreqTimingSignResolution = mode(cell(c).ResolutionTimingSignificativo);
            
            f = figure;
            
            step = intervals/1000 ;
            
            plot((stimulus_start_time/1000):step:(stimulus_end_time(c)/1000),entropy,'b');
            title(cell(c).label);
            ylim([0 1]);

            hold on;

            plot((stimulus_start_time/1000):step:(stimulus_end_time(c)/1000),Qidx,'r');

            hold on;

%             s = 1;
%             for i=(stimulus_start_time/1000):step:(stimulus_end_time(c)/1000)
% 
%                 if (cell(c).significativo(s) == 0)
% 
%                     left = i;
%                     right = i + 0.03;
%                     top = entropy(s);
%                     bottom = 0;
%                     x = [left left right right];
%                     y = [bottom top top bottom];
% 
%                       %fill(x,y, 'b','EdgeColor','none');
% 
%                   else
% 
%                     left = i;
%                     right = i + 0.03;
%                     top = entropy(s);
%                     bottom = 0;
%                     x = [left left right right];
%                     y = [bottom top top bottom];
% 
%                       fill(x,y, 'g','EdgeColor','none');
% 
%                 end
% 
%                 s = s + 1;
% 
%             end    

            
            s = 1;
            for i=(stimulus_start_time/1000):step:(stimulus_end_time(c)/1000)

                if cell(c).Timing(s) == 1
                    
                    if cell(c).significativo(s) == 1
                        
                        plot([i i],[entropy(s) entropy(s)],'r.','markersize',20);
                        
                    else
                        
                        plot([i i],[entropy(s) entropy(s)],'ro');
                        
                    end

                else

                    if cell(c).significativo(s) == 1
                        
                        plot([i i],[entropy(s) entropy(s)],'g.','markersize',20);
                        
                    else
                        
                        plot([i i],[entropy(s) entropy(s)],'go');
                        
                    end

                end

                hold on;

                s = s + 1;

            end

            % p1 = [1 1];
            % p2 = [0 size(Qidx)];
            % plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);
            
            print(f,'-depsc',strcat(cell(c).label,'-resolution-along-train'));
            
            close all;
        
            origem(c).Protocol = [];
            
    end
    
    filepath = 'cells-resolutions-along-train';
    
    save(filepath,'cell','-v7.3');

    toc
    
end

