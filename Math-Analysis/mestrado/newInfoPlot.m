function newInfoPlot

Protocolo(1).metric = load('v3-WOSP-12-08-31-sitio1-E1-full-movie-original.mat');

Protocolo(1).label = 'v3-WOSP-12-08-31-sitio1-E1-full-movie-original';

Protocolo(2).metric = load('v3-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');

Protocolo(2).label = 'v3-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1';

Protocolo(3).metric = load('v3-SP-12-08-31-sitio1-E1-0-10000-full-movie-ae.mat');

Protocolo(3).label = 'v3-SP-12-08-31-sitio1-E1-0-10000-full-movie-ae';

Protocolo(4).metric = load('v311-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');

Protocolo(4).label = 'v311-WOSP-12-12-17-sitio1-E1-full-movie-original';

Protocolo(5).metric = load('v311-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1.mat');

Protocolo(5).label = 'v311-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1';

Protocolo(6).metric = load('v311-SP-12-12-17-sitio1-E1-0-10000-full-movie-ae.mat');

Protocolo(6).label = 'v311-SP-12-12-17-sitio1-E1-0-10000-full-movie-ae';

Protocolo(7).metric = load('v3-WOSP-12-08-29-sitio1-E1-full-movie-original.mat');

Protocolo(7).label = 'v3-WOSP-12-08-29-sitio1-E1-full-movie-original';

Protocolo(8).metric = load('v3-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1.mat');

Protocolo(8).label = 'v3-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1';

Protocolo(9).metric = load('v3-SP-12-08-29-sitio1-E1-0-10000-full-movie-ae.mat');

Protocolo(9).label = 'v3-SP-12-08-29-sitio1-E1-0-10000-full-movie-ae';

Protocolo(10).metric = load('v2-WOSP-12-09-04-sitio1-E1-full-movie-original.mat');

Protocolo(10).label = 'v2-WOSP-12-09-04-sitio1-E1-full-movie-original';

Protocolo(11).metric = load('v2-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1.mat');

Protocolo(11).label = 'v2-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1';

Protocolo(12).metric = load('v2-SP-12-09-04-sitio1-E1-0-10000-full-movie-ae.mat');

Protocolo(12).label = 'v2-SP-12-09-04-sitio1-E1-0-10000-full-movie-ae';


for m=0:1

    for i=1:length(Protocolo)

        % OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        opts.entropy_estimation_method = {'plugin','tpmc','jack'};
        %opts.variance_estimation_method = {'jack'};

        %opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
        opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
        %opts.unoccupied_bins_strategy = 1; % Use all bins

        opts.parallel = 1;
        opts.possible_words = 'unique';

        opts.start_time = 0 / 1000;
        opts.end_time = 10000 / 1000;
        opts.shift_cost = [0 2.^(-4:9)];
        %opts.label_cost = [0 1 2];
        opts.clustering_exponent = -2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        shift_cost = Protocolo(i).metric.metric_analysis.opts.shift_cost;
        X = Protocolo(i).metric.metric_analysis.X;

        opts.entropy_estimation_method = {'plugin'};
        opts.metric_family = m;
        [out_unjk,jk,opts_used] = metric_jack(X,opts);

        P_total = size(jk,1);

        temp_info_jk = zeros(P_total,length(shift_cost));
        
        if m == 0
            
            out = Protocolo(i).metric.metric_analysis.out_T;
            HBias = Protocolo(i).metric.metric_analysis.max_info.HBias_T;
            HBias_std = Protocolo(i).metric.metric_analysis.max_info.HBias_std_T;
            
        elseif m == 1
                
            out = Protocolo(i).metric.metric_analysis.out_I;
            HBias = Protocolo(i).metric.metric_analysis.max_info.HBias_I;
            HBias_std = Protocolo(i).metric.metric_analysis.max_info.HBias_std_I;
            
        end
            
            
        for q_idx=1:length(shift_cost)

            info_tpmc(q_idx) = out(q_idx).table.information(2).value;
            info_jack(q_idx) = out(q_idx).table.information(3).value;
          
            info_unjk(q_idx)= out_unjk(q_idx).table.information.value;

            for p=1:P_total

                temp_info_jk(p,q_idx) = jk(p,q_idx).table.information.value;

            end

        end

        info_jk_sem = sqrt((P_total-1)*var(temp_info_jk,1,1));

        f = figure;
        if m == 0, color = 'b'; else color = 'r'; end
        plot(1:length(shift_cost),info_tpmc,color);
        hold on;
        errorbar(1:length(shift_cost),info_unjk,2*info_jk_sem,color);
        hold on;

        for k=1:length(shift_cost)

            if ( info_tpmc(k) >  ( HBias(k) + 2*HBias_std(k) ) )

                if m == 0, color = 'r.'; else color = 'b.'; end
                plot([k k],[info_tpmc(k) info_tpmc(k)],color,'markersize',20);

            else

                if m == 0, color = 'ro'; else color = 'bo'; end
                plot([k k],[info_tpmc(k) info_tpmc(k)],color);

            end

            hold on;

        end

        set(gca,'xtick',1:length(shift_cost));
        set(gca,'xticklabel',shift_cost);
        set(gca,'xlim',[1 length(shift_cost)]);

        if m == 0, entropy = 'timing'; else entropy = 'interval'; end
        print(f,'-depsc',strcat(Protocolo(i).label,'-info-entropy-',entropy));   

        close all;

    end

end

end

