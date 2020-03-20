function allFractalFFForBack

tic

registro = importdata('memoryBackwardProtocols.txt');

getfractal(registro);

    function getfractal(registro)

        nConditions = 3;

        start_time = 500;

        end_time = 9500;

        bins = 1000;

        for r=1:length(registro)
       
            name = char(registro(r));
            
            if r < 27
                
                protocol = name(2:end-7);
                
            else
                
                protocol = name(2:end-9);
                
            end
            
            Forward_across_bins = load(char(strcat(protocol,'-FF-across_bins-Forward.mat')));
            
            fractalFanoForBack(Forward_across_bins,'Forward','across_bins');
            
            close all;
            
            clear Forward_acroos_bins;
            
            Forward_across_trials = load(char(strcat(protocol,'-FF-across_trials-Forward.mat')));
            
            fractalFanoForBack(Forward_across_trials,'Forward','across_trials');
            
            close all;
            
            clear Forward_across_trials;
                        
            Backward_across_bins = load(char(strcat(protocol,'-FF-across_bins-Backward.mat')));
            
            fractalFanoForBack(Backward_across_bins,'Backward','across_bins');
            
            close all;
            
            clear Backward_across_bins;
            
            Backward_across_trials = load(char(strcat(protocol,'-FF-across_trials-Backward.mat')));
            
            fractalFanoForBack(Backward_across_trials,'Backward','across_trials');
            
            close all;
            
            clear Backward_across_trials;
            
            AE_across_bins = load(char(strcat(protocol,'-FF-across_bins-AE.mat')));
            
            fractalFanoForBack(AE_across_bins,'_AE','across_bins');
            
            close all;
            
            clear AE_acroos_bins;
            
            AE_across_trials = load(char(strcat(protocol,'-FF-across_trials-AE.mat')));
            
            fractalFanoForBack(AE_across_trials,'_AE','across_trials');
            
            close all;
            
            clear AE_across_trials;
            
        end

    end
        
end

